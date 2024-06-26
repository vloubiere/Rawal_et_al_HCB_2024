This repository contains all custom R scripts supporting the conclusion of the study entitled "Sustained activation of the Polycomb PRC1 complex induces DNA repair defects and genomic instability in epigenetic tumors".

Authors: Chetan C. Rawal*, Vincent Loubiere*, Nadedja L. Butova, Julietter Garcia, Victoria Parreno, Chiara Merigliano, Anne-Marie Martinez# & Giacomo Cavalli#, and Irene Chiolo#. (\* Contributed equally; \# Corresponding authors)

System requirements:

All the custom scripts generated for this study were written in R (version 4.2.0) using the R studio IDE (<https://www.R-project.org/>). All custom functions were wrapped into a R package that was made publicly available at <https://github.com/vloubiere/vlfunction/tree/nature_v2_revised>. No special hardware should be required.

Installation guide:

R and RStudio can be downloaded at <https://posit.co/download/rstudio-desktop/>. Installation time is approximately 20min. The R package containing custome function can be installed using the following command: devtools::install_github("<https://github.com/vloubiere/vlfunction/tree/nature_v2_revised>"). Installation time is approximately 20min.

Instructions for use:

The "main.R" file lists all the scripts that are needed to generate the figures of the article. Script names should be self-explanatory, with additional comments linking them to the corresponding figure panel when relevant.

All raw and processed NGS data were retrieved from GEO (GSE222193) and eventually pre-processed using the scripts available at <https://github.com/vloubiere/Parreno_Loubiere_2023.git>. Further details regarding the pre-processing of the data can be found <https://doi.org/10.1038/s41586-024-07328-w>.

For any reasonable further request, please contact [chiolo\@ucsc.edu](mailto:chiolo@ucsc.edu).
