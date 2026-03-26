from dataclasses import dataclass, asdict
from typing import Optional


@dataclass
class VariantAnnotation:
    """Standardized representation for both mtDNA and nucDNA annotations."""

    # Core parsed attributes (Populated during initial parsing)
    ann_id: str  # e.g., 'm.8573G>A' or ClinVar AlleleID
    locus: str  # e.g., 'MT-ND5' or 'SDHA'
    nc_change: str  # e.g., 'm.8573G>A' or 'c.327G>C'
    aa_change: str  # e.g., 'L109F'
    is_synonymous: bool  # True if synonymous, False if missense
    disease: str  # Phenotype or disease string
    genome: str  # 'mtDNA' or 'nucDNA'
    clinical_status: str  # Original status string from the database

    # Database-specific metadata (Populated during initial parsing)
    ref_nt: str
    alt_nt: str
    genomic_pos: int
    clinvar_stars: Optional[int] = None
    clinvar_review_status: Optional[str] = None

    # Downstream classification attributes (Populated later by the Classifier)
    tier: Optional[str] = None
    pathogenic_score: Optional[float] = None

    def to_dict(self):
        """Serializes the dataclass into a dictionary for JSON/CSV export."""
        return asdict(self)
