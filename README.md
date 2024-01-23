# HSC cis inhibition Notch1 Jag1
This repository includes scripts required for the scRNAseq and RNAseq data analysis included in Thambyrajah et al. (2024).
All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

### Main figures folder
Scripts required to reproduce main paper figures, specifically:
- Figure 2: scRNAseq data visualization
- Figure 3: scRNAseq data visualization
- Figure 6: RNAseq data visualization

In order to execute the code, source (processed) data must be downloaded from GEO repository:

- Accession ID GSE230792 for scRNAseq data. Download Cell annotation (cell metadata) and expression matrix (norm counts).
- Accession ID GSE230790 for RNAseq data. Download the expression matrix (normalized counts)

### RNAseq data analysis folder
Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Scripts from 1 to 5.
- Downstream analysis: to conduct differential expression analysis and functional analysis. Scripts 6 and 7.

To conduct data preprocessing, original FASTQ files are required. Please check GEO repository ID GSE230790. All required
scripts were executed using Singularity images (v3.8.3) per each required tool.

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded.

NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 102, mm10).
