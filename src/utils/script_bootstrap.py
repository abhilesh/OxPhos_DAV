from __future__ import annotations

import sys
from pathlib import Path


def bootstrap_src_path(script_file: str) -> None:
    sys.path.insert(0, str(Path(script_file).resolve().parents[1]))
