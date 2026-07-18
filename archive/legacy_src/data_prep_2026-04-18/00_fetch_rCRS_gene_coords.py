#!/usr/bin/env python3
import runpy
from pathlib import Path

if __name__ == "__main__":
    target = Path(__file__).resolve().parents[1] / "data_download" / "00g_fetch_mtdna_gene_coords.py"
    runpy.run_path(str(target), run_name="__main__")
