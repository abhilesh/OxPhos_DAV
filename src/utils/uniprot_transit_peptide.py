"""
src/utils/uniprot_transit_peptide.py

Fetch transit peptide annotations from UniProt REST API for human mitochondrial
proteins.  Used to calibrate mature-protein offset entries in the structural
anchor exception registry.

The transit peptide end position (1-based, inclusive) from UniProt gives the
exact number of N-terminal residues cleaved upon mitochondrial import.  This
offset is used to reconcile RefSeq/TOGA coordinates (pre-protein) with PDB
chain coordinates (mature protein).
"""

import json
import time
import urllib.request
import urllib.error
from typing import Optional

_BASE = "https://rest.uniprot.org/uniprotkb"
_FIELDS = "accession,gene_names,sequence,ft_transit,ft_chain,ft_signal"

# In-process cache: gene → result dict (or None)
_CACHE: dict[str, Optional[dict]] = {}


def _get_json(url: str, retries: int = 3, delay: float = 1.5) -> Optional[dict]:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            if attempt < retries - 1:
                time.sleep(delay)
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def fetch_transit_annotation(
    gene: str,
    organism_id: int = 9606,
    use_cache: bool = True,
) -> Optional[dict]:
    """
    Return transit peptide annotation for a human gene from Swiss-Prot (reviewed).

    Returns a dict with:
      accession       : UniProt accession (e.g. "P49821")
      gene_name       : primary HGNC gene name as reported by UniProt
      sequence        : full canonical pre-protein sequence (one-letter)
      transit_end     : 1-based position of last transit-peptide residue (None if absent)
      mature_start    : transit_end + 1 (= 1 if no transit peptide)
      offset          : transit_end (= 0 if no transit peptide); this is the number
                        of AA to subtract from a pre-protein position to get the
                        mature-protein position
      evidence        : evidence tag string from UniProt (e.g. "ECO:0000269|PubMed")
      feature_type    : "Transit peptide" | "Signal peptide" | None

    Returns None if no reviewed entry is found or the API is unreachable.
    """
    cache_key = f"{gene}:{organism_id}"
    if use_cache and cache_key in _CACHE:
        return _CACHE[cache_key]

    url = (
        f"{_BASE}/search"
        f"?query=gene_exact%3A{gene}+AND+organism_id%3A{organism_id}+AND+reviewed%3Atrue"
        f"&format=json&fields={_FIELDS}&size=1"
    )
    data = _get_json(url)
    if not data or not data.get("results"):
        _CACHE[cache_key] = None
        return None

    entry = data["results"][0]
    accession = entry.get("primaryAccession", "")

    # Primary gene name from UniProt
    gene_name = ""
    for gn_block in entry.get("genes", []):
        gn_val = gn_block.get("geneName", {})
        if isinstance(gn_val, dict):
            for item in gn_val.get("values", []):
                v = item.get("value", "")
                if v:
                    gene_name = v
                    break
        if gene_name:
            break

    seq = entry.get("sequence", {}).get("value", "")

    transit_end: Optional[int] = None
    evidence_str: Optional[str] = None
    feature_type: Optional[str] = None

    # Prefer "Transit peptide"; fall back to "Signal peptide" (secreted subunits)
    for wanted_type in ("Transit peptide", "Signal peptide"):
        for feat in entry.get("features", []):
            if feat.get("type") != wanted_type:
                continue
            loc = feat.get("location", {})
            end_val = loc.get("end", {}).get("value")
            if end_val is not None:
                transit_end = int(end_val)
                feature_type = wanted_type
                evs = feat.get("evidences", [])
                if evs:
                    ev = evs[0]
                    code = ev.get("evidenceCode", "")
                    src_val = ev.get("source")
                    if isinstance(src_val, dict):
                        source = src_val.get("name", "")
                        src_id = src_val.get("id", "")
                    elif isinstance(src_val, str):
                        source, src_id = src_val, ""
                    else:
                        source, src_id = "", ""
                    parts = [p for p in [code, source, src_id] if p]
                    evidence_str = "|".join(parts) if parts else "manual"
                else:
                    evidence_str = "manual_no_code"
                break
        if transit_end is not None:
            break

    result = {
        "accession": accession,
        "gene_name": gene_name,
        "sequence": seq,
        "transit_end": transit_end,
        "mature_start": (transit_end + 1) if transit_end else 1,
        "offset": transit_end if transit_end else 0,
        "evidence": evidence_str,
        "feature_type": feature_type,
    }
    _CACHE[cache_key] = result
    return result
