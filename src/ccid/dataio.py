from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_matrix(path: str | Path, index_col: int | str = 0) -> pd.DataFrame:
    path = Path(path)
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".parquet"):
        frame = pd.read_parquet(path)
        if isinstance(index_col, str) and index_col in frame.columns:
            frame = frame.set_index(index_col)
        elif isinstance(index_col, int) and frame.index.name is None:
            frame = frame.set_index(frame.columns[index_col])
    else:
        sep = "," if ".csv" in suffixes else "\t"
        header = pd.read_csv(path, sep=sep, nrows=0)
        if isinstance(index_col, int):
            index_name = str(header.columns[index_col])
        else:
            index_name = str(index_col)
        dtype_map = {
            str(column): np.float32
            for column in header.columns
            if str(column) != index_name
        }
        frame = pd.read_csv(
            path,
            sep=sep,
            index_col=index_col,
            low_memory=False,
            dtype=dtype_map,
        )
    frame.index = frame.index.astype(str)
    if frame.index.duplicated().any():
        numeric = frame.apply(pd.to_numeric, errors="coerce")
        frame = numeric.groupby(level=0).mean()
    return frame


def write_metadata_template(sample_ids: list[str], output_path: str | Path) -> None:
    pd.DataFrame(
        {
            "sample_id": sample_ids,
            "donor_id": "",
            "day": "",
            "condition": "",
            "cell_type": "",
            "batch": "",
        }
    ).to_csv(output_path, index=False)
