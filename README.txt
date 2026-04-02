## Code and Data for Information-Content-Informed Kendall-tau Correlation Methodology: Interpreting Missing Values in Metabolomics as Potentially Useful Information

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18625643.svg)](https://doi.org/10.5281/zenodo.18625643)

This repository contains the code and partial data for the manuscript:

**Information-Content-Informed Kendall-tau Correlation Methodology: Interpreting Missing Values in Metabolomics as Potentially Useful Information**, R.M. Flight, P.S. Bhatt, and H.N.B. Moseley.


## Manuscripts

* [Main Manuscript - HTML](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTauMetabolomics/ici_kt_manuscript.html)
* [Main Manuscript - Word](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTauMetabolomics/ici_kt_manuscript_mdpi.docx)
* [Supplemental Materials - HTML](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTauMetabolomics/supplemental_materials.html)
* [Supplemental Materials - Word](https://moseleybioinformaticslab.github.io/manuscript.ICIKendallTauMetabolomics/supplemental_materials.docx)

## Grant Support

This work was supported by the following grants:

* NSF 2020026 (PI Moseley);
* NIH 1R03LM014928-01 (PI Moseley);
* NSF ACI1626364 (Griffioen, Moseley);
* P30 CA177558 (PI Evers) via the Markey Cancer Center Biostatistics and Bioinformatics Shared Resource Facility (MCC BB-SRF);
* P20 GM121327 (PD St. Clair); 
* P42 ES007380 (PI Pennell) via the Data Management and Analysis Core (DMAC).

## License

The contents of this work are licensed under a CC-BY license.
If you use any content, you must give attribution to this original work.

## R Packages Needed

The repository used R 4.4.1, and {renv} 1.0.10
For all of the other packages needed, see the file `renv.lock`.

To setup to be able to rerun everything here, you can clone the repo from [github](https://github.com/MoseleyBioinformaticsLab/manuscript.ICIKendallTauMetabolomics) or download it from [Zenodo](https://doi.org/10.5281/zenodo.18625643):

```
# clone from github
git clone https://github.com/MoseleyBioinformaticsLab/manuscript.ICIKendallTauMetabolomics.git
```

```
# download from zenodo
wget 'https://zenodo.org/records/19112314/files/MoseleyBioinformaticsLab/manuscript.ICIKendallTauMetabolomics-v_0.3.zip?download=1' --output-document=manuscript.ICIKendallTauMetabolomics.zip
unzip manuscript.ICIKendallTauMetabolomics.zip
```

```
cd manuscript.ICIKendallTauMetabolomics
```

Start an R session, and make sure `renv` is installed.

```r
# make sure renv is installed
install.packages("renv")
# restore the packages
renv::restore()
```

The {ICIKendallTau} package on GitHub you want to use is the 1.2.16 release.
This is available on zenodo as well.

```r
renv::install("moseleybioinformaticslab/ICIKendallTau@v_1.2.16")
```

Or:

```
wget 'https://zenodo.org/records/18675151/files/MoseleyBioinformaticsLab/ICIKendallTau-v_1.2.16.zip?download=1' --output-document=ICIKendallTau-v_1.2.16.zip
unzip ICIKendallTau-v_1.2.16.zip
cd ICIKendallTau-v_1.2.16
```

Then start an R session from within the directory, and install it.

```r
remotes::install_local()
```

## Obtaining targets Cache and Supporting Datasets

To rerun the analyses, you will need at least the Metabolomics Workbench (MW) datasets and the compound annotations, all hosted on Zenodo.
We also include the instructions to download the original `targets` cache so you can avoid recomputing everything if you like.

```
# required
mkdir mwtab
cd mwtab
wget https://zenodo.org/records/19115655/files/mwbench.tgz?download=1 --output-document=mwbench.tgz
tar -xzf mwbench.tgz
cd ..

# required
mkdir predicted_annotations
cd predicted_annotations
wget https://zenodo.org/records/19115655/files/predicted_annotations.tgz?download=1 --output-document=predicted_annotations.tgz
tar -xzf predicted_annotations.tgz
cd ..

# optional, but suggested
# Note this is 26GB, it will take some time to download!
cd manuscript.ICIKendallTauMetabolomics
wget https://zenodo.org/records/19115655/files/icikt_metabolomics_targets.tgz?download=1 --output-document=icikt_metabolomics_targets.tgz
tar -xzf icikt_metabolomics_targets.tgz
cd ..
```

From within the manuscript project, you will want to link to the MW datasets and the annotations.

```
cd manuscript.ICIKendallTauMetabolomics
ln -s mwtab path/to/mwtab
ln -s predicted_annotations path/to/predicted_annotations/annotations
```

Now you can start R, and hopefully nothing is outdated.

```r
targets::tar_source(c("R", "./packages.R"))
tar_outdated()
```

## Rerunning

This project consists of 91,669 targets for the analysis.
This is because I really wasn't thinking when I designed it, and didn't think to group the MW datasets into batches to be processed.
It made my life easier, but that is a **lot** of things, even though 61,637 are `NULL` objects that still have to be checked whenever you do anything involving the workflow beyond loading things.

### Manuscript and Supplemental Materials

This is why the manuscript and supplemental materials **are not** actual targets of the workflow.
If you just want to check them, then you can do:

```r
rmarkdown::render("docs/supplemental_materials.Rmd", output_format = "rmarkdown::html_document", knit_root_dir = getwd())
rmarkdown::render("docs/ici_kt_manuscript.Rmd", output_format = "rmarkdown::html_document", knit_root_dir = getwd())
```

### Everything Else

You might want to uncomment at the top of the `_targets.R` file to use some multiprocessing, and watch your RAM useage.

```r
targets::tar_make()
```
