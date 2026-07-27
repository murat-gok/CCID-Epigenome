from __future__ import annotations

from itertools import combinations

import numpy as np

from .objectives import evaluate_mask, repair_mask, scalarized_score
from .schema import OptimizationContext, SolutionRecord


def constrained_dominates(a: SolutionRecord, b: SolutionRecord) -> bool:
    if a.feasible and not b.feasible:
        return True
    if not a.feasible and b.feasible:
        return False
    if not a.feasible and not b.feasible:
        return a.violations < b.violations
    return bool(np.all(a.objectives <= b.objectives) and np.any(a.objectives < b.objectives))


def nondominated_front(records: list[SolutionRecord]) -> list[SolutionRecord]:
    keep: list[SolutionRecord] = []
    for i, record in enumerate(records):
        dominated = False
        for j, other in enumerate(records):
            if i != j and constrained_dominates(other, record):
                dominated = True
                break
        if not dominated:
            keep.append(record)
    return _deduplicate_records(keep)


def fast_nondominated_sort(records: list[SolutionRecord]) -> list[list[int]]:
    domination_sets: list[list[int]] = [[] for _ in records]
    domination_counts = np.zeros(len(records), dtype=int)
    fronts: list[list[int]] = [[]]

    for p in range(len(records)):
        for q in range(p + 1, len(records)):
            if constrained_dominates(records[p], records[q]):
                domination_sets[p].append(q)
                domination_counts[q] += 1
            elif constrained_dominates(records[q], records[p]):
                domination_sets[q].append(p)
                domination_counts[p] += 1

    fronts[0] = np.flatnonzero(domination_counts == 0).tolist()
    rank = 0
    while rank < len(fronts) and fronts[rank]:
        next_front: list[int] = []
        for p in fronts[rank]:
            for q in domination_sets[p]:
                domination_counts[q] -= 1
                if domination_counts[q] == 0:
                    next_front.append(q)
        rank += 1
        if next_front:
            fronts.append(next_front)
    return fronts


def crowding_distance(records: list[SolutionRecord], indices: list[int]) -> dict[int, float]:
    if len(indices) <= 2:
        return {idx: float("inf") for idx in indices}
    values = np.vstack([records[idx].objectives for idx in indices])
    distances = np.zeros(len(indices), dtype=float)
    for objective in range(values.shape[1]):
        order = np.argsort(values[:, objective])
        distances[order[0]] = np.inf
        distances[order[-1]] = np.inf
        minimum = values[order[0], objective]
        maximum = values[order[-1], objective]
        if maximum <= minimum:
            continue
        for position in range(1, len(indices) - 1):
            prev_value = values[order[position - 1], objective]
            next_value = values[order[position + 1], objective]
            distances[order[position]] += (next_value - prev_value) / (maximum - minimum)
    return {indices[i]: float(distances[i]) for i in range(len(indices))}


def _deduplicate_records(records: list[SolutionRecord]) -> list[SolutionRecord]:
    unique: dict[bytes, SolutionRecord] = {}
    for record in records:
        key = np.packbits(record.mask.astype(np.uint8)).tobytes()
        current = unique.get(key)
        if current is None or scalarized_score(record.objectives, record.violations) < scalarized_score(
            current.objectives, current.violations
        ):
            unique[key] = record
    return list(unique.values())


def _make_record(
    mask: np.ndarray,
    context: OptimizationContext,
    method: str,
    seed: int,
) -> SolutionRecord:
    objectives, feasible, violations = evaluate_mask(mask, context)
    return SolutionRecord(
        mask=mask.astype(bool),
        objectives=objectives,
        feasible=feasible,
        violations=violations,
        method=method,
        seed=seed,
    )


def random_search(
    context: OptimizationContext,
    n_samples: int = 2500,
    seed: int = 0,
) -> list[SolutionRecord]:
    rng = np.random.default_rng(seed)
    records: list[SolutionRecord] = []
    n_candidates = len(context.candidate_table)
    for _ in range(n_samples):
        k = int(rng.integers(context.min_edits, context.max_edits + 1))
        mask = np.zeros(n_candidates, dtype=bool)
        if k > 0:
            mask[rng.choice(n_candidates, size=k, replace=False)] = True
        mask = repair_mask(mask, context, rng)
        records.append(_make_record(mask, context, "random", seed))
    return nondominated_front(records)


def greedy_search(context: OptimizationContext, seed: int = 0) -> list[SolutionRecord]:
    rng = np.random.default_rng(seed)
    n_candidates = len(context.candidate_table)
    current = np.zeros(n_candidates, dtype=bool)
    records = [_make_record(current, context, "greedy", seed)]

    for _ in range(context.max_edits):
        candidates: list[SolutionRecord] = []
        for idx in np.flatnonzero(~current):
            trial = current.copy()
            trial[idx] = True
            trial = repair_mask(trial, context, rng)
            candidates.append(_make_record(trial, context, "greedy", seed))
        if not candidates:
            break
        best = min(candidates, key=lambda r: scalarized_score(r.objectives, r.violations))
        if np.array_equal(best.mask, current):
            break
        current = best.mask.copy()
        records.extend(candidates)
    return nondominated_front(records)


def exact_enumeration(
    context: OptimizationContext,
    candidate_indices: np.ndarray | None = None,
    seed: int = 0,
) -> list[SolutionRecord]:
    n_candidates = len(context.candidate_table)
    if candidate_indices is None:
        candidate_indices = np.arange(n_candidates)
    candidate_indices = np.asarray(candidate_indices, dtype=int)
    if len(candidate_indices) > 18:
        raise ValueError("Exact enumeration is intentionally limited to <=18 candidates")

    records: list[SolutionRecord] = []
    for k in range(context.min_edits, context.max_edits + 1):
        for combo in combinations(candidate_indices.tolist(), k):
            mask = np.zeros(n_candidates, dtype=bool)
            if combo:
                mask[list(combo)] = True
            records.append(_make_record(mask, context, "exact", seed))
    return nondominated_front(records)


def nsga2(
    context: OptimizationContext,
    population_size: int = 96,
    generations: int = 90,
    crossover_rate: float = 0.90,
    mutation_rate: float | None = None,
    seed: int = 0,
) -> list[SolutionRecord]:
    rng = np.random.default_rng(seed)
    n_candidates = len(context.candidate_table)
    if mutation_rate is None:
        mutation_rate = min(0.20, 2.0 / n_candidates)

    population: list[SolutionRecord] = []
    for _ in range(population_size):
        mask = rng.random(n_candidates) < (context.max_edits / max(n_candidates, 1))
        mask = repair_mask(mask, context, rng)
        population.append(_make_record(mask, context, "nsga2", seed))

    for _ in range(generations):
        fronts = fast_nondominated_sort(population)
        rank_map: dict[int, int] = {}
        crowd_map: dict[int, float] = {}
        for rank, front in enumerate(fronts):
            for idx in front:
                rank_map[idx] = rank
            crowd_map.update(crowding_distance(population, front))

        def tournament() -> SolutionRecord:
            i, j = rng.integers(0, len(population), size=2)
            if rank_map[i] < rank_map[j]:
                return population[i]
            if rank_map[j] < rank_map[i]:
                return population[j]
            if crowd_map.get(i, 0.0) > crowd_map.get(j, 0.0):
                return population[i]
            return population[j]

        offspring: list[SolutionRecord] = []
        while len(offspring) < population_size:
            parent_a = tournament().mask
            parent_b = tournament().mask
            child_a = parent_a.copy()
            child_b = parent_b.copy()
            if rng.random() < crossover_rate:
                selector = rng.random(n_candidates) < 0.5
                child_a = np.where(selector, parent_a, parent_b)
                child_b = np.where(selector, parent_b, parent_a)
            child_a ^= rng.random(n_candidates) < mutation_rate
            child_b ^= rng.random(n_candidates) < mutation_rate
            child_a = repair_mask(child_a, context, rng)
            child_b = repair_mask(child_b, context, rng)
            offspring.append(_make_record(child_a, context, "nsga2", seed))
            if len(offspring) < population_size:
                offspring.append(_make_record(child_b, context, "nsga2", seed))

        combined = _deduplicate_records(population + offspring)
        fronts = fast_nondominated_sort(combined)
        next_population: list[SolutionRecord] = []
        for front in fronts:
            if len(next_population) + len(front) <= population_size:
                next_population.extend(combined[idx] for idx in front)
            else:
                distances = crowding_distance(combined, front)
                ordered = sorted(front, key=lambda idx: distances[idx], reverse=True)
                slots = population_size - len(next_population)
                next_population.extend(combined[idx] for idx in ordered[:slots])
                break
        population = next_population

    return nondominated_front(population)


def select_compromise(front: list[SolutionRecord]) -> SolutionRecord:
    feasible = [
        record
        for record in front
        if record.feasible and record.objectives[2] >= 1.0
    ]
    pool = feasible if feasible else front
    matrix = np.vstack([record.objectives for record in pool])
    minimum = matrix.min(axis=0)
    span = np.maximum(matrix.max(axis=0) - minimum, 1e-12)
    normalized = (matrix - minimum) / span
    distance = np.sqrt(np.sum(np.square(normalized), axis=1))
    return pool[int(np.argmin(distance))]


def select_target_priority(front: list[SolutionRecord]) -> SolutionRecord:
    """Prespecified decision rule for the primary biological endpoint.

    Feasibility is mandatory; target matching receives dominant weight, while
    safety, uncertainty and burden break near-ties.
    """
    feasible = [
        record
        for record in front
        if record.feasible and record.objectives[2] >= 1.0
    ]
    pool = feasible if feasible else front
    matrix = np.vstack([record.objectives for record in pool])
    minimum = matrix.min(axis=0)
    span = np.maximum(matrix.max(axis=0) - minimum, 1e-12)
    normalized = (matrix - minimum) / span
    score = (
        0.72 * normalized[:, 0]
        + 0.16 * normalized[:, 1]
        + 0.05 * normalized[:, 2]
        + 0.07 * normalized[:, 3]
    )
    return pool[int(np.argmin(score))]


def binary_differential_evolution(
    context: OptimizationContext,
    population_size: int = 96,
    generations: int = 90,
    crossover_rate: float = 0.80,
    mutation_rate: float | None = None,
    seed: int = 0,
) -> list[SolutionRecord]:
    """Multiweight binary differential-evolution baseline.

    The binary differential operator uses XOR differences. Selection is based
    on randomly sampled Tchebycheff decompositions, with feasibility taking
    precedence. The returned archive is the constrained nondominated set.
    """

    rng = np.random.default_rng(seed)
    n_candidates = len(context.candidate_table)
    n_objectives = 4
    if population_size < 4:
        raise ValueError("population_size must be at least 4")
    if mutation_rate is None:
        mutation_rate = min(0.15, 1.0 / max(n_candidates, 1))

    weights = rng.dirichlet(np.ones(n_objectives), size=population_size)
    population: list[SolutionRecord] = []
    for _ in range(population_size):
        mask = rng.random(n_candidates) < (context.max_edits / max(n_candidates, 1))
        mask = repair_mask(mask, context, rng)
        population.append(_make_record(mask, context, "binary_de", seed))
    archive = population.copy()

    def decomposition(record: SolutionRecord, weight: np.ndarray, ideal: np.ndarray) -> float:
        if not record.feasible:
            return 1e6 + 1e4 * record.violations
        return float(np.max(np.maximum(weight, 1e-6) * np.abs(record.objectives - ideal)))

    for _ in range(generations):
        feasible_values = [record.objectives for record in population if record.feasible]
        ideal = (
            np.min(np.vstack(feasible_values), axis=0)
            if feasible_values
            else np.min(np.vstack([record.objectives for record in population]), axis=0)
        )
        next_population = population.copy()
        for i, target in enumerate(population):
            pool = np.delete(np.arange(population_size), i)
            a_idx, b_idx, c_idx = rng.choice(pool, size=3, replace=False)
            a = population[a_idx].mask
            b = population[b_idx].mask
            c = population[c_idx].mask
            donor = a ^ (b ^ c)
            trial = target.mask.copy()
            cross = rng.random(n_candidates) < crossover_rate
            if n_candidates:
                cross[int(rng.integers(0, n_candidates))] = True
            trial[cross] = donor[cross]
            trial ^= rng.random(n_candidates) < mutation_rate
            trial = repair_mask(trial, context, rng)
            child = _make_record(trial, context, "binary_de", seed)

            if constrained_dominates(child, target):
                next_population[i] = child
            elif not constrained_dominates(target, child):
                if decomposition(child, weights[i], ideal) <= decomposition(target, weights[i], ideal):
                    next_population[i] = child
        population = next_population
        archive = nondominated_front(archive + population)
    return nondominated_front(archive)


def moead(
    context: OptimizationContext,
    population_size: int = 96,
    generations: int = 90,
    neighborhood_size: int = 15,
    crossover_rate: float = 0.90,
    mutation_rate: float | None = None,
    seed: int = 0,
) -> list[SolutionRecord]:
    """Compact MOEA/D baseline for the four-objective binary protocol problem."""

    rng = np.random.default_rng(seed)
    n_candidates = len(context.candidate_table)
    n_objectives = 4
    if population_size < 4:
        raise ValueError("population_size must be at least 4")
    neighborhood_size = int(np.clip(neighborhood_size, 2, population_size))
    if mutation_rate is None:
        mutation_rate = min(0.20, 2.0 / max(n_candidates, 1))

    weights = rng.dirichlet(np.ones(n_objectives), size=population_size)
    distances = np.linalg.norm(weights[:, None, :] - weights[None, :, :], axis=2)
    neighborhoods = np.argsort(distances, axis=1)[:, :neighborhood_size]

    population: list[SolutionRecord] = []
    for _ in range(population_size):
        mask = rng.random(n_candidates) < (context.max_edits / max(n_candidates, 1))
        mask = repair_mask(mask, context, rng)
        population.append(_make_record(mask, context, "moead", seed))
    archive = nondominated_front(population)

    def tchebycheff(record: SolutionRecord, weight: np.ndarray, ideal: np.ndarray) -> float:
        if not record.feasible:
            return 1e6 + 1e4 * record.violations
        return float(np.max(np.maximum(weight, 1e-6) * np.abs(record.objectives - ideal)))

    feasible_values = [record.objectives for record in population if record.feasible]
    ideal = (
        np.min(np.vstack(feasible_values), axis=0)
        if feasible_values
        else np.min(np.vstack([record.objectives for record in population]), axis=0)
    )

    for _ in range(generations):
        for i in rng.permutation(population_size):
            neighborhood = neighborhoods[i]
            p1, p2 = rng.choice(neighborhood, size=2, replace=False)
            parent_a = population[p1].mask
            parent_b = population[p2].mask
            child_mask = parent_a.copy()
            selector = rng.random(n_candidates) < crossover_rate
            child_mask[selector] = parent_b[selector]
            child_mask ^= rng.random(n_candidates) < mutation_rate
            child_mask = repair_mask(child_mask, context, rng)
            child = _make_record(child_mask, context, "moead", seed)
            if child.feasible:
                ideal = np.minimum(ideal, child.objectives)

            for j in neighborhood:
                incumbent = population[j]
                if constrained_dominates(child, incumbent) or (
                    not constrained_dominates(incumbent, child)
                    and tchebycheff(child, weights[j], ideal)
                    <= tchebycheff(incumbent, weights[j], ideal)
                ):
                    population[j] = child
            archive = nondominated_front(archive + [child])
    return nondominated_front(archive)
