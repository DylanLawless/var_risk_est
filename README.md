# Quantifying the Genetics of Disease Inheritance for Bayesian Application

VarRiskEst is an integrative framework that quantifies the probability of observing disease-associated variants across the genome by combining large-scale genomic annotation data with classical Hardy–Weinberg equilibrium (HWE) calculations. By integrating data from gnomAD, ClinVar, dbNSFP, and curated gene panels from PanelAppRex, our method estimates observation probabilities for each single nucleotide variant (SNV) and aggregates these probabilities by gene and variant classification. These estimates serve as informative priors for Bayesian models of variant and disease risk, offering enhanced precision in genetic diagnostics.

## Data Sources

- **PanelAppRex**: Aggregates disease gene panel data from sources such as GE PanelApp, ClinVar, and UniProt to facilitate natural language-based searches for clinical and research applications \[Lawless et al., 2025\].
- **ClinVar**: A public archive of human genetic variation that provides detailed variant classifications and supporting evidence \[Landrum et al., 2018\].
- **dbNSFP**: A comprehensive database for the functional prediction and annotation of non-synonymous and splicing-site SNVs in human protein-coding genes \[Liu et al., 2020\].
- **gnomAD**: Provides reference allele frequencies from diverse human populations.

## Overview

Our study focused on a disease gene panel (ID 398) associated with primary immunodeficiency or monogenic inflammatory bowel disease. We analyzed 54,814 ClinVar variant classifications across 557 genes to estimate the probability of disease observation. The framework accounts for complexities in different inheritance modes (autosomal recessive, autosomal dominant, and X-linked) by incorporating both homozygous/compound heterozygous states and single pathogenic alleles into our HWE-based calculations.

The resulting observation probabilities not only refine classical risk estimates but also serve as robust priors in a Bayesian framework. This approach enriches clinical understanding and supports bioinformatic analyses, ultimately leading to more precise genetic diagnostics.

## Installation

The project requires R (version X.X or later) along with several R packages. To install the required packages, run:

```r
install.packages(c("dplyr", "readr", "patchwork", "stringr", "knitr", "kableExtra"))
```

Clone the repository:

```bash
git clone https://github.com/DylanLawless/var_risk_est.git
cd var_risk_est
```

## Usage

The main computational pipeline is implemented in R scripts within the repository. Key scripts include:

- **main_pipeline.R**: Aggregates variant data, computes allele frequencies, and estimates observation probabilities.
- **validate_nfkb1.R**: Validates autosomal dominant estimates using *NFKB1*.
- **validate_cftr.R**: Validates autosomal recessive estimates using *CFTR* (p.Arg117His).
- **PanelAppRex_import.R**: Imports and processes PanelAppRex data.

To run the main analysis, start R and execute:

```r
source("main_pipeline.R")
```

Results (summary tables and figures) will be saved in the `../images/` and `../output/` directories.

## Output Files

- **Data Files**: Final processed data are saved in both TSV and RDS formats, named as follows:
  - `VarRiskEst_<PanelAppRex_ID>_gene_variants.tsv` and `VarRiskEst_<PanelAppRex_ID>_gene_variants.Rds`
  - `VarRiskEst_<PanelAppRex_ID>_gene_tally.tsv` and `VarRiskEst_<PanelAppRex_ID>_gene_tally.Rds`
- **Figures**: Validation study figures and probability distribution plots are stored in the `../images/` directory.

## Citation

Please consider citing the following references when using VarRiskEst:

```
@article {Lawless2025VarRiskEst,
	author = {Lawless, Dylan},
	title = {Quantitative prior probabilities for disease-causing variants reveal the top genetic contributors in inborn errors of immunity},
	elocation-id = {2025.03.25.25324607},
	year = {2025},
	doi = {10.1101/2025.03.25.25324607},
	publisher = {Cold Spring Harbor Laboratory Press}
}

```

Please consider citing the related works used to develop this tool:

- Lawless, D. et al. (2025). PanelAppRex aggregates disease gene panels and facilitates sophisticated search. [GitHub repository](https://github.com/DylanLawless/PanelAppRex).
- Landrum, M. J. et al. (2018). ClinVar: updates to support classifications of both germline and somatic variants. *Nucleic Acids Res.* doi:10.1093/nar/gkae1090.
- Liu, X. et al. (2020). dbNSFP v4.0: A comprehensive database of functional predictions and annotations for human nonsynonymous and splice-site SNVs. *Nucleic Acids Res.*
- Martin, A. R. et al. (2019). PanelApp: a tool for the diagnosis of rare diseases. *Nat. Genet.* [Include full citation details].


## License

VarRiskEst is available under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

We acknowledge GE PanelApp for providing public access to curated disease gene panels, as well as ClinVar and UniProt for their indispensable data resources. We also thank our collaborators for their feedback and support.

