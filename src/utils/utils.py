import sys
from pathlib import Path


def get_cached_path(target_path: Path, date_str: str) -> Path:
    """
    Checks if a previously downloaded version of the file exists by replacing
    the date string with a wildcard '*' and finding the most recent match.
    """
    pattern = target_path.name.replace(date_str, "*")
    matches = list(target_path.parent.glob(pattern))

    if matches:
        return max(matches, key=lambda p: p.stat().st_mtime)

    return target_path


def make_read_only(file_path: Path):
    """Locks the file to prevent accidental overwrite by curation scripts."""
    if file_path.exists():
        # 0o444 gives Read permission to User, Group, and Others (no write access)
        file_path.chmod(0o444)
        print(f"    Locked {file_path.name} as Read-Only.")


def get_latest(data_dir: Path, pattern: str) -> Path:
    """
    Recursively searches a directory for the most recently modified file
    matching a specific glob pattern.
    """
    try:
        return max(data_dir.rglob(pattern), key=lambda f: f.stat().st_mtime)
    except ValueError:
        print(
            f"Error: File missing matching '{pattern}'. Please run 00b_download_annotation_data.py"
        )
        sys.exit(1)
