"""Antibody Analysis Tools Subpackage.

This subpackage provides various tools for analyzing antibody sequences, including:
- Sequence homology searches (BLAST)
- Physicochemical property calculation (ProtParam)
- Stability prediction
- Aggregation propensity prediction
- Glycosylation site prediction
- Structure prediction (via SWISS-MODEL)
- B-cell epitope prediction (via IEDB)
- Epitope conservancy analysis
- Developability assessment based on sequence similarity to known antibodies.
"""

# Import primary functions and classes from submodules to make them easily accessible

from .blast_analyzer import perform_blast
from .protparam_analyzer import analyze_with_protparam
from .stability_predictor import predict_stability
from .aggregation_predictor import predict_aggregation
from .glycosylation_analyzer import predict_glycosylation
from .structure_predictor import predict_structure
from .epitope_predictor import predict_bcell_epitopes, get_epitope_sequences_from_prediction, VALID_IEDB_BCELL_METHODS
from .conservancy import EpitopeConservancyAnalyzer, predict_conservancy
from .developability import AntibodyComparisonTool, predict_developability
from .predict_affinity import predict_binding_affinity
from .immunogenicity import predict_immunogenicity

# Define __all__ for explicit public API
__all__ = [
    # Functions
    "perform_blast",
    "analyze_with_protparam",
    "predict_stability",
    "predict_aggregation",
    "predict_glycosylation",
    "predict_structure",
    "predict_bcell_epitopes",
    "get_epitope_sequences_from_prediction", # Helper, but potentially useful
    "VALID_IEDB_BCELL_METHODS", # 
    "predict_conservancy",
    "predict_developability",
    "predict_binding_affinity",
    "predict_immunogenicity",
    
    # Classes
    "EpitopeConservancyAnalyzer",
    "AntibodyComparisonTool",
]