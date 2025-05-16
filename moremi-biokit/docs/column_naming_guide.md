# Protein Ranking Column Naming Guide

This document explains the column naming convention used in the protein ranking output CSV files.

## Column Types

The ranking output includes several types of columns for each metric category:

1. **Raw Values**: The actual measurements or values calculated for each metric
2. **Normalized Scores**: Values transformed to a 0-1 scale based on optimal ranges
3. **Category Scores**: Combined/averaged scores for a category of metrics
4. **Weighted Scores**: Category scores multiplied by their importance weight

## Naming Conventions

### 1. Raw Values

Raw metric values are named with a descriptive term and `_value` suffix for clarity:

- `gravy_value`: The GRAVY value calculated from ProtParam
- `gmqe_value`: The GMQE score from structure prediction
- `immunogenicity_value`: The raw immunogenicity score
- `developability_value`: The raw developability score

Some special metrics retain their original names for clarity:

- `melting_temperature`: The melting temperature in Celsius
- `solubility`: Categorical value (e.g., "Soluble")
- `aggregation_propensity`: Categorical value (e.g., "Low", "Medium", "High")
- `aggregation_regions_count`: Count of aggregation-prone regions
- `n_glyc_sites_count`: Count of N-glycosylation sites

### 2. Normalized Scores

Normalized scores (0-1 scale) begin with `norm_` prefix:

- `norm_gravy_score`: Normalized GRAVY score (0-1)
- `norm_solubility_score`: Normalized solubility score (0-1)
- `norm_agg_propensity_score`: Normalized aggregation propensity score (0-1)
- `norm_agg_regions_score`: Normalized aggregation regions count score (0-1)

### 3. Category Scores

Category scores have the category name with `_score` suffix:

- `protparam_score`: Overall score for the ProtParam category
- `immunogenicity_score`: Overall score for the Immunogenicity category
- `stability_score`: Overall score for the Stability category
- `aggregation_score`: Overall score for the Aggregation category
- `glycosylation_score`: Overall score for the Glycosylation category
- `structure_score`: Overall score for the Structure category
- `binding_affinity_score`: Overall score for the Binding Affinity category
- `epitope_score`: Overall score for the Epitope category
- `conservancy_score`: Overall score for the Conservancy category
- `developability_score`: Overall score for the Developability category

### 4. Weighted Scores

Weighted scores begin with `weighted_` prefix:

- `weighted_protparam_score`: ProtParam score × weight
- `weighted_immunogenicity_score`: Immunogenicity score × weight
- `weighted_stability_score`: Stability score × weight
- ...and so on for all categories

## Additional Important Columns

- `sequence`: The protein sequence
- `antigen`: The antigen sequence
- `antigen_pdb_chain_id`: PDB chain ID of the antigen (e.g., "4R19_A" or "predicted_structure")
- `molecular_formula`: Molecular formula of the protein
- `molecular_weight`: Molecular weight of the protein
- `total_score`: Final combined score (0-1 scale) - higher is better
- `warnings`: Number of warnings raised during processing
- `warning_details`: Details of warnings
- `missing_metrics`: Number of metrics that couldn't be calculated

## Example

For the ProtParam category:

- `gravy_value`: -0.352 (actual calculated GRAVY value)
- `norm_gravy_score`: 0.852 (normalized to 0-1 scale based on optimal range)
- `protparam_score`: 0.926 (average of normalized GRAVY and solubility scores)
- `weighted_protparam_score`: 0.0926 (category score × weight of 0.10)

## Interpreting Scores

- **Raw values** provide actual measurements for detailed analysis
- **Normalized scores** show how good each metric is on a 0-1 scale
- **Category scores** combine related metrics into a single 0-1 score
- **Weighted scores** account for the importance of each category
- **Total score** is the sum of all weighted category scores
