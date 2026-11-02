# Survival-Associated Transcriptomic Signatures in Glioblastoma

Code, filtered data and supplementary tables for my BMedSci (Honours) Bioinformatics research project at the University of Cape Town (2026).

**Project title:** Survival-Associated Transcriptomic Signatures in Glioblastoma: Differential Expression Analysis and Alignment with Established Molecular Subtypes

**Author:** Njabulo Dlomo
**Supervisor:** Dr Hocine Bendou
**Co-supervisor:** Dr Michelle Livesey

## Overview

This project used bulk RNA-seq data from 122 TCGA glioblastoma patients, all with a recorded death, to look for genes associated with survival. Patients were grouped into short survival (SS, < 365 days), intermediate survival (IS, 365–730 days) and long survival (LS, > 730 days). The analysis has three parts:

1. **Differential expression and enrichment.** DESeq2 contrasts between the three survival groups, followed by GO and KEGG over-representation analysis and Hallmark GSEA.
2. **Survival analysis.** Univariate Cox regression, log-rank tests, LASSO-Cox feature selection with stability selection, a six-gene risk score, and evaluation with repeated cross-validation and time-dependent ROC.
3. **Molecular subtypes.** Overlap of the survival-associated genes with the Verhaak (2010) and Crisman (2016) subtype signatures, the association between subtype and survival group, and an elastic-net classifier testing whether genes outside the Verhaak list can classify the classical, mesenchymal and proneural subtypes.

## Repository contents

| Folder / file | Contents |
|---|---|
| `GBM_Honours_Pipeline.R` | Full analysis pipeline, from data filtering to the final figures and tables |
| `Filtered Files/` | Filtered expression and phenotype files used as input to the survival and subtype analyses |
| `Crisman's & Verhaak's Supplementary Files/` | Verhaak ClaNC centroid file and Crisman 48-gene panel (see Data sources) |
| `Supplementary Tables/` | Supplementary Tables S1–S8 from the thesis |

## Data sources

- **Expression and clinical data:** GDC TCGA-GBM, downloaded from the [UCSC Xena Browser](https://xenabrowser.net/) (Goldman et al., 2020). Gene-level STAR counts and open-access phenotype data only.
- **Molecular subtype labels:** retrieved with the TCGAbiolinks R package (Colaprico et al., 2016) from the annotation table of Ceccarelli et al. (2016).
- **Verhaak subtype signatures:** the ClaNC 840-gene centroid file from the [GDC publication page](https://gdc.cancer.gov/about-data/publications/gbm_exp) (Verhaak et al., 2010). The four subtype gene sets were also retrieved from MSigDB using msigdbr.
- **Crisman 48-gene panel:** extracted from Table 2 of Crisman et al. (2016).
- **Gene annotation:** protein-coding genes defined with biomaRt (Ensembl release 116, GRCh38.p14).

## How to run

1. Download or clone this repository.
2. Open `GBM_Honours_Pipeline.R` in RStudio.
3. Set the working directory to the repository folder, or update the file paths at the top of the script.
4. Install the packages listed below if needed.
5. Run the script from top to bottom. The sections follow the order of the thesis Methods.

Random seeds are set in the script (42 for GSEA and the main LASSO-Cox fit, 1–50 for the LASSO-Cox stability fits, 2026 and 1001–1010 for the subtype classifier) so that the results can be reproduced.

## Software

Analyses were run in R 4.5.2 using RStudio.

| Package | Version | Used for |
|---|---|---|
| biomaRt | 2.66.2 | Protein-coding gene annotation |
| DESeq2 | 1.50.2 | Normalisation, VST and differential expression |
| ggplot2 | 4.0.2 | Figures |
| ggVennDiagram | 1.5.7 | Venn diagram |
| EnhancedVolcano | 1.28.2 | Volcano plots |
| clusterProfiler | 4.18.4 | GO, KEGG and GSEA |
| org.Hs.eg.db | 3.22.0 | Gene annotation for enrichment |
| msigdbr | 26.1.0 | Hallmark and Verhaak gene sets |
| survival | 3.8-3 | Cox regression and log-rank tests |
| survminer | 0.5.2 | Kaplan–Meier plots |
| glmnet | 5.0 | LASSO-Cox and elastic-net models |
| timeROC | 0.4.1 | Time-dependent ROC |
| TCGAbiolinks | 2.38.0 | Molecular subtype annotations |
| caret | 7.0.1 | Stratified cross-validation folds |
| HGNChelper | – | Updating Verhaak gene symbols |

## Notes and limitations

- All patients in the cohort had died, so there is no censoring, and patients still alive at last follow-up are not represented.
- The six-gene risk score was developed and evaluated in the same TCGA cohort and has not been validated in an independent cohort. The results should be read as hypothesis-generating.

## Citation

If you use this code, please cite the project:

Dlomo, N. (2026). *Survival-associated transcriptomic signatures in glioblastoma: Differential expression analysis and alignment with established molecular subtypes* [BMedSci (Honours) research project, University of Cape Town].

## Licence

The code is released under the MIT licence. The Verhaak and Crisman gene lists remain the work of their original authors; please cite the original publications when using them:

- Verhaak, R. G. W., et al. (2010). Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1. *Cancer Cell*, 17(1), 98–110. https://doi.org/10.1016/j.ccr.2009.12.020
- Crisman, T. J., et al. (2016). Identification of an efficient gene expression panel for glioblastoma classification. *PLoS ONE*, 11. https://doi.org/10.1371/journal.pone.0164649
