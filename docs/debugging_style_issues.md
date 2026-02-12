---
title: 'Information-Content-Informed Kendall-tau Correlation Methodology: Interpreting Missing Values in Metabolomics as Potentially Useful Information'
author:
  - Robert M Flight:
      institute: [markey, biochem, rcsirm]
  - Praneeth S Bhatt:
      institute: compeng
  - Hunter NB Moseley:
      email: hunter.moseley@uky.edu
      correspondence: true
      institute: [markey, biochem, rcsirm, ibi, tox]
institute:
  - markey: Markey Cancer Center, University of Kentucky, Lexington, KY 40536, USA
  - biochem: Department of Molecular & Cellular Biochemistry, University of Kentucky, Lexington, KY 40536, USA
  - rcsirm: Resource Center for Stable Isotope Resolved Metabolomics, University of Kentucky, Lexington, KY 40536, USA
  - compeng: Department of Electrical and Computer Engineering, University of Kentucky, Lexington, KY 40506
  - ibi: Institute for Biomedical Informatics, University of Kentucky, Lexington, KY 40536, USA
  - tox: Department of Toxicology and Cancer Biology, University of Kentucky, Lexington, KY 40536, USA
output: 
  rmarkdown::word_document:
    keep_md: true
    reference_docx: metabolites-template-alt.dot
    pandoc_args:
      - '--lua-filter=scholarly-metadata.lua'
      - '--lua-filter=author-info-blocks.lua'
  rmarkdown::html_document:
    self_contained: yes
    pandoc_args:
      - '--lua-filter=scholarly-metadata.lua'
      - '--lua-filter=author-info-blocks.lua'
date: '2026-02-12 10:54:13.349671'
bibliography: '/home/rmflight/Documents/manuscripts/in_progress/rmflight_icikt_metabolomics/docs/icikt_manuscript.json'
csl: plos-computational-biology.csl
editor_options: 
  chunk_output_type: console
---









::: {custom-style="MDPI_1.7_abstract"}
**Abstract**

**Background:** 
Almost all correlation measures currently available are unable to directly handle missing values.
Typically, missing values are either ignored completely by removing them or are imputed and used in the calculation of the correlation coefficient.
In either case, the correlation value will be impacted based on a perspective that the missing data represents no useful information.
However, missing values occur in real data sets for a variety of reasons. 
In metabolomics data sets, a major reason for missing values is that a specific measurable phenomenon falls below the detection limits of the analytical instrumentation (left-censored values).
These missing data are not missing at random, but represent potentially useful information by virtue of their "missingness" at one end of the data distribution.
**Methods:** 
To include this information due to left-censorship missingness, we propose the information-content-informed Kendall-tau (ICI-Kt) methodology.
We develop a statistical test showing that most missing values are the result of left-censorship, and show how left-censored missing values can be included within the definition of the Kendall-tau correlation coefficient, and how that inclusion leads to an interpretation of information being added to the correlation.
We also implement calculations for additional measures of theoretical maxima and pairwise completeness that add further layers of information interpretation in the methodology.
**Results:** 
Using both simulated and over 700 data sets from metabolomics repositories, we demonstrate that the ICI-Kt methodology allows for the inclusion of left-censored missing data values as interpretable information, enabling both improved determination of outlier samples and improved feature-feature network construction.
**Conclusions:** 
We provide explicitly parallel implementations in both R and Python that allow fast calculations of all the variables used when applying the ICI-Kt methodology on large numbers of samples.
The ICI-Kt methods are available as an R package and Python module on GitHub at https://github.com/moseleyBioinformaticsLab/ICIKendallTau and https://github.com/moseleyBioinformaticsLab/icikt, respectively.
:::

::: {custom-style="MDPI_1.8_keywords"}
**Keywords:** Metabolomics; Correlation; Missingness; Left-Censored
:::

::: {custom-style="MDPI_2.1_heading1"}
**1. Introduction**
:::

::: {custom-style="MDPI_3.1_text"}
Correlation as a measure of the relatedness or similarity of two or more sets of data has a long history, with the mathematical technique being used (and abused) in various scientific fields since its introduction [@pearson_notes_1920; @rodgers_13_1988].
More recently, correlation calculations have become a cornerstone statistical method in the analysis and integration of varied omics' datasets, especially the big five omics: genomics, transcriptomics, proteomics, metabolomics, and epigenomics [@gu_complexheatmap_2016].
Correlation between biomolecular features (nucleotide variants, RNA transcripts, proteins, metabolites, DNA methylation, etc.) may be used to evaluate the relationship strength between pairs of the features as well as to detect and derive correlative structures between groups of features [@fukushima_integratedomics_2009].
Moreover, feature-feature correlations can be used to evaluate a dataset based on expected biochemical correlations, for example higher feature-feature correlations within lipid categories versus between lipid categories [@mitchellUntargetedLipidomicsNonSmall2021].
Correlation is a foundational method for generating biomolecular feature-feature interaction networks, like those provided by STRING [@szklarczyk_string_2017], Genemania [@franz_genemania_2018], and WCGNA [@langfelder_wgcna_2008].
Feature-feature correlation may also be used to inform which features are used for imputation of missing values [@faquih_missingvalueworkflow_2020].

Often, the first step in omics level analyses is to examine the sample-sample (dis)similarities in various ways using exploratory data analysis or EDA.
This can include the examination of decomposition by principal components analysis (PCA), sample-sample pairwise distances, or sample-sample pairwise correlations to highlight biological and batch groups [@loveRNASeqWorkflowGenelevel2016; @lawRNAseqAnalysisEasy2018; @chenReadsGenesPathways2016], double check the appropriateness of planned analyses [@flight_timecourseexploration_2010], and check if any samples should be removed prior to statistical analysis (outlier detection and removal) [@gierlinski_statisticalmodels_2015].
Outlier detection, in particular, is often required for successful omics data analysis, as any misstep during the experimentation, sample collection, sample preparation, or analytical measurement of individual samples can inject high error and/or variance into the resulting data set [@loveRNASeqWorkflowGenelevel2016; @lawRNAseqAnalysisEasy2018; @chenReadsGenesPathways2016; @moseleyErrorAnalysisPropagation2013; @gierlinski_statisticalmodels_2015].


All analytical methods, and in particular the analytical methods used in omics where many analytes are being measured simultaneously, suffer from missing measurements.
Some analytes will be missing at random because of spurious issues with either the instrument,  the particular sample, or sample preparation, but a larger number of missing measurements are left-censored due to analytes being below the effective detection limit of the instrument and the given specific sample preparation procedures utilized, as shown in Figure 1.
Some analytical instruments are purposely designed to floor measurements when they occur below a certain signal to noise ratio threshold.
Also, imputation of missing measurements in omics samples is an active area of research, which we will not comprehensively cover here beyond to say that it is worthwhile and very necessary in many instances.
Imputation methods rely on very similar analytical detection limits between analytical samples. 
When this condition does not hold, imputation methods have reduced performance and lower interpretive value. 
For analytical techniques requiring complex sample handling and detection, the variability in the analytical detection level can be quite high.
However, when it comes to calculating correlation, there are very few methods that explicitly account for left-censored missing data that we know of.
In many cases, missing values are either ignored or imputed to zero (or another value) and then included in the correlation calculation.
The two most common approaches for ignoring (i.e. dropping) values is to only use those measurements that are common across all samples (complete) or that are common between two samples being compared (pairwise complete).
Both dropping or imputing missing values are likely to cause the calculated sample-sample correlation values to deviate from the real sample-sample correlation values, especially with respect to specific data interpretation perspectives.
:::

::: {custom-style="MDPI_5.2_figure"}
![](/home/rmflight/Documents/manuscripts/in_progress/rmflight_icikt_metabolomics/docs/the_problem.png){width=400 height=403}
:::

::: {custom-style="MDPI_5.1_figure_caption"}
**Figure 1.** 
Graphical description of the left-censored data problem.
An example density plot of the analyte concentrations for a single experimental sample is shown as a solid black curve. 
The true analyte concentration range covers the full range of the density distribution, with the minimum on the left (red vertical line), and the maximum on the right (yellow vertical line).
Below certain concentrations, shown by the red line, the instrument returns either missing values (NA), zeros, or some other floored values, resulting in a left-censored distribution.
Above certain concentrations, highlighted by the yellow line, typically the values returned will be identical (or flagged depending on the instrument).
Which analytes will have concentrations below the red detection limit line may vary from sample to sample due to the overall sample composition, as well as natural variances (non-systematic error) within each experimental sample.
:::

::: {custom-style="MDPI_3.1_text"}
Assuming that a majority of missing values are not missing at random, but rather result from left-censored distributions due to the analyte being below the effective detection limit (see Figure 1), we propose that these missing values do in fact encode useful information that can be incorporated into correlation calculations.

To create a correlation measure that is capable of working with missing values, we would not be interested in creating a completely new correlation metric from scratch, but modifying an existing one.
Of the three commonly used correlation measures, Pearson, Spearman, and Kendall-tau, Spearman and Kendall-tau seem most appropriate for modification as they solely use ranks in the calculation of their coefficients.
Modifying Pearson would either involve imputing new values, or finding a way to calculate the covariances **with** missingness included.
While Spearman uses ranks, many of the modifications for handling identical ranks and ties do not seem amenable to working with missing values.
In contrast, Kendall-tau's use of concordant and discordant pair counts seems most amenable to the creation of new definitions that incorporate missingness while still working within the original definition of the correlation coefficient, as shown in the Implementation section below.

In this work, we propose new definitions of concordant and discordant rank pairs that include missing values, as well as methods for incorporating missing values into the number of tied values for the equivalent of the modified Kendall-tau coefficient, the information-content-informed Kendall-tau (ICI-Kt) method.
The implementation of the basic calculation of ICI-Kt involves the replacement of missing values with a value lower than the observed values (technically simple imputation), with subsequent calculation of the Kendall-tau-b statistic, as a majority of missing values are the result of left-censorship, they still provide an interpretation from an information-content perspective, which we demonstrate with the equations below.
We also developed a statistical test for determining if the cause for missingness is likely left-censorship. With this statistical test, we demonstrate that left-censorship is likely the cause of many missing values across a large number of metabolomics datasets from the Metabolomics Workbench (MW).
We examine the effect of missing values on various collections of simulated and real data-sets, comparing the ICI-Kt methodology with other simpler methods of handling the missing values, namely removing them or imputing them to zero.
Given the detrimental effects of including outlier samples in differential analyses, we also evaluate the ability of the ICI-Kt methodology to capture sample-sample pairwise similarities and the determination of outlier samples prior to differential analyses.
We were also curious about the utility of the ICI-Kt methodology in creating feature-feature networks with large amounts of missing values, so we evaluate the partitioning of networks by reactome pathways after network creation using different correlation measures.

All of the code used for this manuscript is available on zenodo [@flightManuscriptICIKendallTau2024].
:::

::: {custom-style="MDPI_2.1_heading1"}
**2. Materials and Methods**
:::

::: {custom-style="MDPI_2.2_heading2"}
*2.1 Additional Definitions of Concordant and Discordant Pairs to Include Missingness*
:::

::: {custom-style="MDPI_3.1_text"}
In the simplest form, the Kendall-tau correlation can be defined as:
:::

::: {custom-style="MDPI_3.9_equation"}
$$\tau_a = \frac{ n_{concordant} - n_{discordant} }{n_{concordant} + n_{discordant}}$$ [1]{custom-style="MDPI_3.a_equation_number"}
:::

::: {custom-style="MDPI_3.1_text"}
where $n_{concordant}$ is the the number of concordant pairs and $n_{discordant}$ is the number of discordant pairs. 
In this case, a pair are any two x-y points, $x_i, y_i$ and $x_j, y_j$, with $i \neq j$, composed from two joint random variables X and Y, where $x_i$ represents the *ith value* in X and $y_i$ represents the *ith value* in Y.
In an omics context, X and Y can represent feature vectors for two experimental samples or two specific features across a set of samples.
:::
::: {custom-style="MDPI_3.5_text_before_list"}
A concordant pair has the following classical definition:
:::
::: {custom-style="MDPI_3.8_bullet_alt"}
$x_i \gt x_j$ and $y_i \gt y_j$;

$x_i \lt x_j$ and $y_i \lt y_j$;
:::

::: {custom-style="MDPI_3.5_text_before_list"}
A discordant pair has the following classical definition [@kendall_newmeasure_1938]:
:::

::: {custom-style="MDPI_3.8_bullet_alt"}
$x_i \gt x_j$ and $y_i \lt y_j$

$x_i \lt x_j$ and $y_i \gt y_j$
:::

We can expand the concordant and discordant pair definitions to include missing values (e.g. `NA` in R). 
Here $!x$ indicates $x=\text{NA}$ and $!y$ indicates $y=\text{NA}$, and $\&$ is synonymous with "and".
The information-content-informed concordant pair definitions are then:

  * $x_i > x_j$ and $y_i > y_j$
  * $x_i < x_j$ and $y_i < y_j$
  * $x_i > x_j$ and $y_i \& !y_j$
  * $x_i < x_j$ and $!y_i \& y_j$
  * $x_i \& !x_j$ and $y_i > y_j$
  * $!x_i \& x_j$ and $y_i < y_j$ 
  * $x_i \& !x_j$ and $y_i \& !y_j$ 
  * $!x_i \& x_j$ and $!y_i \& y_j$ 
  
The information content informed discordant pair definitions are then:

  * $x_i > x_j$ and $y_i < y_j$ 
  * $x_i < x_j$ and $y_i > y_j$
  * $x_i > x_j$ and $!y_i \& y_j$ 
  * $x_i < x_j$ and $y_i \& !y_j$
  * $x_i \& !x_j$ and $y_i < y_j$
  * $!x_i \& x_j$ and $y_i > y_j$
  * $x_i \& !x_j$ and $!y_i \& y_j$
  * $!x_i \& x_j$ and $y_i \& !y_j$
  
