#!/usr/bin/env python3
import runpy
from pathlib import Path

if __name__ == "__main__":
    target = Path(__file__).resolve().parents[1] / "data_download" / "00e_build_taxid_species_map.py"
    runpy.run_path(str(target), run_name="__main__")
