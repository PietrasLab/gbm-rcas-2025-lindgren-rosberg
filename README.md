# GBM RCAS – single-cell and spatial transcriptomics analysis

This repository contains analysis code, metadata, and documentation for the study:

**Coevolution of neoplastic and non-neoplastic reactive astrocyte states converges on mesenchymal-like and injury-response programs during murine glioblastoma progression and post-radiotherapy recurrence**

David Lindgren#, Rebecca Rosberg#, Karolina I. Smolag, Elinn Johansson, Dimitra Manou, Jonas Sjölund, Sebastian Braun, Bengt Phung, Katja Harbst, Pauline Jeannot, Maria Malmberg, Crister Ceberg, Göran B. Jönsson, Kristian Pietras, Alexander Pietras*  
\# Equal contribution  
\* Corresponding author

---

## Abstract

Glioblastomas initially respond to radiotherapy but invariably recur, often within high-dose radiation treatment fields. Although stromal radiation responses are incompletely understood, evidence suggests that the tumor microenvironment becomes tumor-supportive after therapy. Using a genetically engineered glioblastoma mouse model, we profiled healthy brain, primary tumors, and post-radiotherapy recurrences with single-cell and spatial transcriptomics and immunohistochemistry. Across 13 non-neoplastic cell types and 10 tumor cell states, we mapped transcriptional adaptations accompanying progression from healthy brain to primary and recurrent glioblastoma. We identified distinct astrocyte states linked to disease stage, including reactive non-neoplastic astrocytes and tumor cells adopting reactive astrocyte-like phenotypes. Reactive astrocyte-like tumor cells with mesenchymal and injury-response signatures were enriched after radiotherapy and persisted in recurrent tumors. Receptor–ligand interactions between reactive astrocytes and tumor cells included known and putative drivers of aggressiveness. These findings highlight convergent reactive astrocyte programs in astrocytes and tumor cells as potential mediators of glioblastoma radioresistance.

---

## Scope of this repository

This repository contains:

- Analysis scripts (R, Python, bash)
- Metadata tables used for analysis and figure generation
- Gene sets, probe definitions, and reference files
- Documentation describing analysis steps and data availability

**Raw sequencing data, primary processing outputs, and large processed objects are intentionally excluded from version control.**

---

## Processed data availability (Zenodo)

Analysis-ready Seurat v5 objects used in this study are deposited on Zenodo:

**Zenodo record:**  
https://zenodo.org/records/21607378

### Key datasets

- **Main processed single-cell dataset (10x Genomics Flex):**  
  `seurat_flex_filtered_v1.0.rds`  
  https://zenodo.org/records/21607378/files/seurat_flex_filtered_v1.0.rds

- **Processed Visium HD spatial transcriptomics datasets:**  
  - Visium HD Seurat objects (8 µm and 16 µm bins; healthy and tumor)

  Visium v1 objects from the initial submission were removed from this deposit and were not used in the final analysis; they remain available in earlier (superseded) versions of this Zenodo record.

- **Raw (unfiltered) single-cell Seurat object:**  
  `seurat_flex_raw_v1.0.rds` 

These objects were used for all downstream analyses and figure generation in the accompanying manuscript.

---

## Raw sequencing and primary processing data

Raw sequencing data (FASTQ files) and primary processing outputs
(Cell Ranger for 10x Genomics Flex and Space Ranger for Visium HD)
have been submitted to ArrayExpress.

These submissions are currently under curation and/or embargo and are not yet public.
Accession numbers will be added here once they are released.

---

## Reproducibility

Analyses were performed primarily in R (Seurat v5 and associated packages), with
additional processing and validation in Python (Scanpy and related tools).
Scripts are organized by analysis stage within the `scripts/` directory.

---

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).