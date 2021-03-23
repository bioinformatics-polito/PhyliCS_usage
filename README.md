# PhyliCS Usage
This repository contains the data and the jupyter notebooks to reproduce the case studies presented in [PhyliCS](https://github.com/bioinformatics-polito/PhyliCS/tree/master) paper.

Spefifically, 
- [synth.ipynb](https://github.com/bioinformatics-polito/PhyliCS_usage/blob/main/synth.ipynb) implements the code to compute SHscore statistics presented in the section **Experiment 1: Synthetic data**
- [breast_tumor.ipynb](https://github.com/bioinformatics-polito/PhyliCS_usage/blob/main/breast_tumor.ipynb) implements the code to reproduce the analysis on 10x Genomics breast tumor data presented in section **Experiment 2: Tumor Data**, subsection **Spatial subsamples from the same disease site**
- [lung_tumor.ipynb](https://github.com/bioinformatics-polito/PhyliCS_usage/blob/main/lung_tumor.ipynb) implements the code to reproduce the analysis on public lung primary tumors and liver metastasis data presented in section **Experiment 2: Tumor Data**, subsection **Spatial subsamples from the different disease sites**
- [data] contains all data needed to reproduce the analyses

In order to read the notebooks content, you just need to click on them and they will automatically be rendered as HTML files. If you prefer, you may clone the current directory and execute the notebooks locally, after [installing PhyliCS](https://github.com/bioinformatics-polito/PhyliCS#installation-and-setup-instructions) and its dependencies.

If you encounter any problem in rendering the notebooks you may navigate to each `data` subfolder where you can find a README describing the use case (e.g. `data\simulations\README.md`).
