# Protein Ranking Output Improvements

## Overview

We've made significant improvements to the protein ranking output to enhance clarity, add important information, and fix inconsistencies.

## Key Changes

### 1. Improved Column Naming Convention

The new naming convention clearly distinguishes between:

- **Raw values**: The actual measurements (e.g., `gravy_value`)
- **Normalized scores**: Values transformed to 0-1 scale (e.g., `norm_gravy_score`)
- **Category scores**: Combined scores for a category (e.g., `protparam_score`)
- **Weighted scores**: Category scores multiplied by weights (e.g., `weighted_protparam_score`)

### 2. Added Antigen PDB Chain ID

- Added `antigen_pdb_chain_id` to both the `ProteinMetrics` class and the ranking output
- This field shows:
  - The PDB and chain ID combination (e.g., "4R19_A") when using an external PDB
  - "predicted_structure" when using a structure predicted from a sequence

### 3. Fixed Inconsistencies in Score Reporting

- Standardized score field names for consistency across metrics
- Changed confusing prefix `raw_` to clearer `norm_` for normalized scores
- Fixed category score naming to clearly indicate the category

## Benefits

- **Easier interpretation**: Clear distinction between raw values and processed scores
- **Better tracking**: Inclusion of antigen PDB chain ID allows better tracking of which structures were used
- **Improved navigation**: Consistent naming makes it easier to find specific metrics
- **Enhanced documentation**: Added a comprehensive column naming guide

## Documentation

For detailed explanations, see:

- [Column Naming Guide](column_naming_guide.md): Detailed explanation of the naming convention
- [Sample Output](../sample_ranking_output_template.csv): Example CSV showing the new format

## Example

Here's an example of the improved naming for the ProtParam metrics:

**Before**:

```plaintext
gravy_score: -0.33
raw_gravy_score: 0.84
raw_protparam_score: 0.92
weighted_protparam_score: 0.092
```

**After**:

```plaintext
gravy_value: -0.33
norm_gravy_score: 0.84
protparam_score: 0.92
weighted_protparam_score: 0.092
```
