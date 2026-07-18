#!/usr/bin/env python3
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notes" / "running_results.md"
RESDIR = ROOT / "data" / "derived" / "results"
CUR = ROOT / "data" / "derived" / "curated" / "variants_master_curated.parquet"
CLS = ROOT / "data" / "derived" / "classified" / "variants_master_classified.parquet"
MANIFEST = ROOT / "data" / "derived" / "reference" / "download_manifest.json"


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            vals.append("" if pd.isna(v) else str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def main() -> None:
    curdf = pd.read_parquet(CUR)
    cldf = pd.read_parquet(CLS)
    manifest_rows = json.loads(MANIFEST.read_text())
    active = [r for r in manifest_rows if r.get("active_for_pipeline")]

    summary_rows = []
    for genome in ["mtDNA", "nucDNA"]:
        sub = cldf[cldf["genome"] == genome]
        subc = sub[sub["classification_status"] == "classified"]
        sube = sub[sub["classification_status"].isin(["classified", "unresolved"])]
        summary_rows.append(
            {
                "Genome": genome,
                "Total rows": len(sub),
                "Classified": len(subc),
                "Eligible": len(sube),
                "Unresolved": int((sub["classification_status"] == "unresolved").sum()),
                "Skipped by policy": int((sub["classification_status"] == "skipped_by_policy").sum()),
                "AA-level cDAVs": int(subc["is_cdav_amino_acid"].fillna(False).sum()),
                "NT-level cDAVs": int(subc["is_cdav_nucleotide"].fillna(False).sum()),
                "AA cDAV % of classified": round(subc["is_cdav_amino_acid"].fillna(False).mean() * 100, 2)
                if len(subc)
                else None,
                "NT cDAV % of classified": round(subc["is_cdav_nucleotide"].fillna(False).mean() * 100, 2)
                if len(subc)
                else None,
            }
        )
    summary_df = pd.DataFrame(summary_rows)

    complex_df = pd.read_csv(RESDIR / "dav_metrics_by_complex.tsv", sep="\t")
    gene_df = pd.read_csv(RESDIR / "dav_metrics_by_gene.tsv", sep="\t")

    hgnc = ROOT / "data" / "raw" / "reference" / "Canonical_OXPHOS_Subunits_HGNC_2026-03-25.csv"
    rows = []
    with hgnc.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            vals = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, vals)))
    syms = [r["Approved symbol"] for r in rows if r.get("Approved symbol")]
    hgnc_total = len(syms)
    hgnc_mt = sum(s.startswith("MT-") for s in syms)
    hgnc_nuc = sum(not s.startswith("MT-") for s in syms)

    cur_mt = curdf[curdf["genome"] == "mtDNA"]
    cur_nuc = curdf[curdf["genome"] == "nucDNA"]

    fixes = [
        "The pipeline was converted to a metadata-first, filter-late framework: records are retained with explicit eligibility and exclusion metadata rather than being dropped during parsing or classification.",
        "Overlapping mtDNA loci, especially MT-ATP6/MT-ATP8, are duplicated into frame-specific interpretations, so one nucleotide event can yield two independent curated and classified rows.",
        'Gene identity was standardized around interpreted_gene, replacing older logic that collapsed composite loci such as locus.split("/")[0].',
        "The nuclear transcript-position map builder was corrected to parse full multi-line TOGA human protein FASTA sequences. This removed a major source of false POSITION_NOT_IN_ENST calls.",
        "Classification now reads the canonical curated Parquet master table, retains all rows, and writes canonical classified outputs under data/derived/classified/.",
        "Parquet-to-JSON export was stabilized by normalizing NumPy arrays and NaN values before writing compatibility JSON and JSONL outputs.",
        "Coordinate rescue is now ordered more defensibly: deterministic transcript maps are preferred when map quality is sufficient, genomic rescue is used for lower-concordance genes, and anchor fallback remains available only as a guarded last-resort heuristic.",
        "Residual problematic genes and rows are now handled via an explicit exception registry plus exported audit tables instead of hidden one-off code paths.",
    ]

    mitigations = [
        "Anchor fallback: if no transcript map exists, the classifier can search within +/-10 amino-acid positions for the expected wild-type residue in the aligned human sequence before projecting codon coordinates.",
        "Genomic rescue: for lower-concordance genes, genomic coordinates are projected directly into the CDS space of the selected TOGA ENST model.",
        "Mismatch categories are explicit and not conflated with uDAV calls: REF_ALLELE_MISMATCH, TRANSCRIPT_MISMATCH, POSITION_NOT_IN_ENST, GENOMIC_POS_NOT_IN_ENST, ANCHOR_NOT_FOUND, CODON_EXTRACTION_FAILURE, NO_ALIGNMENT, and COORD_PARSE_FAILURE.",
        "Classified-with-warning rows are preserved but flagged when a cDAV/uDAV classification succeeds despite a reference-base disagreement.",
    ]

    active_df = pd.DataFrame(
        [
            {
                "Resource": r["resource_name"],
                "Download date": r["download_date"] or "undated_legacy",
                "Local file": r["local_path"],
                "Source URL": r["source_url"],
                "Validation": r["validation_status"],
            }
            for r in active
        ]
    )

    main_complex = complex_df[
        ["genome", "complex_id", "classified_rows", "aa_cdav", "nt_cdav", "aa_cdav_pct", "nt_cdav_pct"]
    ].copy()
    mt_gene_df = gene_df[
        gene_df["genome"] == "mtDNA"
    ][["gene", "complex_id", "classified_rows", "aa_cdav", "nt_cdav", "aa_cdav_pct", "nt_cdav_pct"]].copy()
    nuc_gene_df = gene_df[
        gene_df["genome"] == "nucDNA"
    ][["gene", "complex_id", "classified_rows", "aa_cdav", "nt_cdav", "aa_cdav_pct", "nt_cdav_pct"]].copy()

    text = []
    text.append("# Running Results")
    text.append("")
    text.append(
        "This document records the important analysis results accumulated so far for the OXPHOS DAV study. It will be updated as additional stages of the pipeline are completed."
    )
    text.append("")
    text.append("## 1. Current scope")
    text.append("")
    text.append("- Current implemented stages: `data_download`, `data_curation`, and `classify`.")
    text.append(
        f"- Current targeted OXPHOS gene set from HGNC: `{hgnc_total}` genes total (`{hgnc_mt}` mitochondrial, `{hgnc_nuc}` nuclear)."
    )
    text.append(f"- Current curated DAV inventory: `{len(curdf)}` curated interpretation rows.")
    text.append("- Current disease-variant sources entering curation: `MITOMAP` for mtDNA and `ClinVar` for nucDNA.")
    text.append("")
    text.append("## 2. Data sources used so far")
    text.append("")
    text.append("The current download and reference layer includes the following active source inputs:")
    text.append("")
    text.append(md_table(active_df))
    text.append("")
    text.append("Additional comparative and reference assets currently present in the repo and used by the pipeline include:")
    text.append("")
    text.append(
        "- TOGA codon-aware and amino-acid alignments for nuclear genes under `data/alignments/toga_hg38_aa/` and `data/alignments/toga_hg38_codon/`."
    )
    text.append(
        "- Existing mtDNA codon-aware and amino-acid alignments under `data/alignments/mtdna_codon/` and `data/alignments/mtdna_aa/`."
    )
    text.append("- Canonical transcript and genomic rescue maps under `data/derived/curated/`.")
    text.append("- Exception registry and focused audit products under `data/derived/reference/` and `data/derived/classified/`.")
    text.append("")
    text.append("## 3. Downloaded genes and variant inventory")
    text.append("")
    text.append(f"- HGNC targeted gene list: `{hgnc_total}` OXPHOS genes total (`{hgnc_mt}` mtDNA-encoded, `{hgnc_nuc}` nucDNA-encoded).")
    text.append(
        f"- Curated mtDNA interpretation rows from MITOMAP: `{len(cur_mt)}` rows spanning `{cur_mt['interpreted_gene'].nunique()}` genes."
    )
    text.append(
        f"- Curated nucDNA interpretation rows from ClinVar: `{len(cur_nuc)}` rows spanning `{cur_nuc['interpreted_gene'].nunique()}` genes."
    )
    text.append(
        f"- Unique source-variant groups: `{cur_mt['source_variant_group_id'].nunique()}` for mtDNA and `{cur_nuc['source_variant_group_id'].nunique()}` for nucDNA."
    )
    text.append(
        "- mtDNA curated rows are interpretation-level rows and include overlap-derived duplications where biologically required, especially for `MT-ATP6/MT-ATP8`."
    )
    text.append("")
    text.append("## 4. Current cDAV classification results")
    text.append("")
    text.append(
        "The current classifier identifies compensated and uncompensated DAVs at both the amino-acid and nucleotide levels for mtDNA and nucDNA."
    )
    text.append("")
    text.append(md_table(summary_df))
    text.append("")
    text.append("Main current cDAV proportions using classified rows as the denominator:")
    text.append("")
    for _, row in summary_df.iterrows():
        text.append(
            f"- `{row['Genome']}`: AA-level cDAVs `{row['AA-level cDAVs']}/{row['Classified']}` = `{row['AA cDAV % of classified']}%`; NT-level cDAVs `{row['NT-level cDAVs']}/{row['Classified']}` = `{row['NT cDAV % of classified']}%`."
        )
    text.append("")
    text.append("At the current state of the pipeline, this corresponds to:")
    text.append("")
    text.append("- `mtDNA`: `167` AA-level cDAVs and `149` NT-level cDAVs among `323` classified rows.")
    text.append("- `nucDNA`: `2046` AA-level cDAVs and `1745` NT-level cDAVs among `6486` classified rows.")
    text.append("")
    text.append("## 5. cDAV metrics by OXPHOS complex")
    text.append("")
    text.append(md_table(main_complex))
    text.append("")
    text.append(
        "A tab-delimited copy of the complex-level breakdown is available at [dav_metrics_by_complex.tsv](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/dav_metrics_by_complex.tsv)."
    )
    text.append("")
    text.append("## 6. Supplementary information")
    text.append("")
    text.append("### 6.1 Important fixes and nuanced issues resolved so far")
    text.append("")
    for item in fixes:
        text.append(f"- {item}")
    text.append("")
    text.append("### 6.2 Important mitigations currently in use")
    text.append("")
    for item in mitigations:
        text.append(f"- {item}")
    text.append("")
    text.append("### 6.3 DAV metrics by gene")
    text.append("")
    text.append(
        "The full gene-level breakdown is also exported as [dav_metrics_by_gene.tsv](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/dav_metrics_by_gene.tsv). The current gene-level supplementary tables are reproduced below."
    )
    text.append("")
    text.append("#### mtDNA genes")
    text.append("")
    text.append(md_table(mt_gene_df))
    text.append("")
    text.append("#### nucDNA genes")
    text.append("")
    text.append(md_table(nuc_gene_df))
    text.append("")
    text.append("### 6.4 Current residual classification caveats")
    text.append("")
    text.append(
        "- Residual unresolved nuclear rows remain concentrated in a small set of genes with transcript-model incompatibility or gene-specific transcript/consequence discordance, especially `NDUFS6`, `NDUFA13`, `NDUFA11`, `UQCRB`, `NDUFS7`, `NDUFV2`, `NDUFB1`, and `NDUFA10`."
    )
    text.append("- These are now tracked explicitly via the exception registry rather than being hidden in implicit script logic.")
    text.append(
        "- Focused review tables are available at `data/derived/classified/exception_candidate_rows.tsv` and `data/derived/classified/exception_candidate_summary.tsv`."
    )
    text.append("")

    OUT.write_text("\n".join(text), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
