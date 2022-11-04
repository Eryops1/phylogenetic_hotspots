# Global hotspots of plant phylodiversity

## Summary
We estimate the global patterns of plant phylogenetic diversity and endemism using species-level distribution data for 330,527 seed plant species from the World Checklist for Vascular Plants [1]. We define phylodiversity hotspots as areas that jointly maximise species-count based indices (i.e. species richness and endemism) and standardised phylogenetic indices (i.e. standardised PD and standardised PE). Analysing global plant phylodiversity and identify hotspots, compare with conservation hotspots.

## Content
This repo contains all code necessary to repeat the study. Since some steps include computation-intense steps (e.g. supplementing the phylogeny with TACT, calculation of standardized phylodiversity per tree), we provide some derived data files to enable skipping the data preparation steps.

All analysis were done in R version 4.2.1 (R Core Team, 2022). 

Upon publication, this repo will be archived on Zenodo.

## Notes on computational requirements
### TACT
Adding missing taxa to the phylogenetic tree with TACT [2] takes about 80 hours and requires 150GB RAM for each run. Due to the stochastic nature of the procedure, we ran TACT 100 times, and calculated the phylogeny-derived estimates as average values from all 100 trees.

### sesPD
Each calculation run takes 4h and requires 64GB RAM. We ran this on a cluster and provide here in addition to the R script the bash scripts for SLURM.

## Literature
[1] Govaerts, R., Nic Lughadha, E., Black, N. et al. The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. Sci Data 8, 215 (2021). https://doi.org/10.1038/s41597-021-00997-6

[2] Jonathan Chang, Daniel L Rabosky, Michael E Alfaro, Estimating Diversification Rates on Incompletely Sampled Phylogenies: Theoretical Concerns and Practical Solutions, Systematic Biology, Volume 69, Issue 3, May 2020, Pages 602–611, https://doi.org/10.1093/sysbio/syz081
<!--# Links

[Manuscript](https://docs.google.com/document/d/10n0I9uDsZKgldacGPLwiH3PRePPw8nT8jgoGek2puzQ/edit#)

[Supplement](https://docs.google.com/document/d/1zH4zpbRLBgCsSY9FiQyHNe2FeAxJo6vJO0IdqfL6Mac/edit#)

Notebook for daily work, ideas etc: [google docs](https://docs.google.com/document/d/1xpm09o4z9haBdFzmqCq2ZjAE9qxPqobLEedJZhql3cM/edit#)
-->
