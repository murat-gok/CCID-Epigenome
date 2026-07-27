# Upgrade to CCID-Epigenome v0.9.0.post1

## Install

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.9.0.post1.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected:

```text
64 passed
0.9.0.post1
```

## Run the core orthogonal validation

```bash
bash scripts/run_orthogonal_validation_v090.sh
```

## Run the optional 479-sample lifespan analysis

```bash
RUN_EXTENDED=1 bash scripts/run_orthogonal_validation_v090.sh
```

## Main result

```bash
column -s, -t \
  < results/real/orthogonal_validation_v090/integrated/orthogonally_validated_candidate_panel_v090.csv \
  | less -S
```

## Important scope

v0.9.0.post1 strengthens external validity and biological interpretation using independent public data. It does not replace targeted epigenome-editing experiments and must not be described as wet-lab or causal validation.
