from __future__ import annotations

import gzip
import json
import ssl
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Callable


SSL_CTX = ssl.create_default_context()
SSL_CTX.check_hostname = False
SSL_CTX.verify_mode = ssl.CERT_NONE


def fetch_bytes(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, context=SSL_CTX, timeout=timeout) as response:
        return response.read()


def fetch_stream_to_file(url: str, dest: Path, timeout: int = 60) -> int:
    req = urllib.request.Request(url)
    total = 0
    with urllib.request.urlopen(req, context=SSL_CTX, timeout=timeout) as response, open(dest, "wb") as handle:
        while True:
            chunk = response.read(1024 * 64)
            if not chunk:
                break
            handle.write(chunk)
            total += len(chunk)
    return total


def validate_gzip_tsv(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        next(handle)
    return "ok"


def validate_tsv(path: Path) -> str:
    with open(path, encoding="utf-8", errors="ignore") as handle:
        next(handle)
    return "ok"


def validate_json(path: Path) -> str:
    with open(path, encoding="utf-8") as handle:
        json.load(handle)
    return "ok"


def validate_zip(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        archive.testzip()
    return "ok"


def safe_validate(path: Path, validator: Callable[[Path], str]) -> str:
    try:
        return validator(path)
    except Exception as exc:
        return f"invalid:{exc}"


def remote_probe(url: str, timeout: int = 20) -> dict:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, context=SSL_CTX, timeout=timeout) as response:
            return {
                "reachable": True,
                "status": getattr(response, "status", 200),
                "content_type": response.headers.get("Content-Type", ""),
                "content_length": response.headers.get("Content-Length", ""),
            }
    except urllib.error.HTTPError as exc:
        if exc.code == 405:
            try:
                req = urllib.request.Request(url)
                with urllib.request.urlopen(req, context=SSL_CTX, timeout=timeout) as response:
                    probe = response.read(256)
                    return {
                        "reachable": True,
                        "status": getattr(response, "status", 200),
                        "content_type": response.headers.get("Content-Type", ""),
                        "content_length": response.headers.get("Content-Length", ""),
                        "probe_bytes": len(probe),
                    }
            except urllib.error.HTTPError as nested_exc:
                return {
                    "reachable": False,
                    "status": nested_exc.code,
                    "content_type": nested_exc.headers.get("Content-Type", ""),
                    "content_length": nested_exc.headers.get("Content-Length", ""),
                    "error": str(nested_exc),
                }
            except Exception as nested_exc:
                return {"reachable": False, "error": str(nested_exc)}
        return {
            "reachable": False,
            "status": exc.code,
            "content_type": exc.headers.get("Content-Type", ""),
            "content_length": exc.headers.get("Content-Length", ""),
            "error": str(exc),
        }
    except Exception:
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, context=SSL_CTX, timeout=timeout) as response:
                probe = response.read(256)
                return {
                    "reachable": True,
                    "status": getattr(response, "status", 200),
                    "content_type": response.headers.get("Content-Type", ""),
                    "content_length": response.headers.get("Content-Length", ""),
                    "probe_bytes": len(probe),
                }
        except urllib.error.HTTPError as exc:
            return {
                "reachable": False,
                "status": exc.code,
                "content_type": exc.headers.get("Content-Type", ""),
                "content_length": exc.headers.get("Content-Length", ""),
                "error": str(exc),
            }
        except Exception as exc:
            return {"reachable": False, "error": str(exc)}
