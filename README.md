[![DOI](https://zenodo.org/badge/525328601.svg)](https://zenodo.org/badge/latestdoi/525328601)

# Vides et. al 2022
Companion code for Vides et. al 2022. [[1]](#1)

This repository contains the code used to analyze the TIRF microscopy images in the Rab10-dependent LRRK2 recruitment in vitro assay.

The following steps were performed:
- LRRK2 recruitment movies were processed with the single molecule tracking and analysis framework TrackIT [[2]](#2). 
- Tracks of fluorescent intensity over time of individual LRRK2 spots were exported to csv files using a custom matlab script [tracks_export.m](tracks_export.m)
- Tracks were analyzed in R as reported in the [rab10_tirf.Rmd](Rmd/rab10_tirf.md) notebook.

## References
<a id="1">[1]</a> 
Edmundo G Vides, Ayan Adhikari, Claire Y Chiang, Pawel Lis, Elena Purlyte, Charles Limouse, Justin L Shumate, Elena Spinola-Lasso, Herschel S Dhekne, Dario R Alessi, Suzanne R Pfeffer (2022).
A feed-forward pathway drives LRRK2 kinase membrane recruitment and activation. eLife 11:e79771.
https://doi.org/10.7554/eLife.79771

<a id="2">[2]</a> 
Kuhn, T., Hettich, J., Davtyan, R. et al. 
Single molecule tracking and analysis framework including theory-predicted parameter settings. Sci Rep 11, 9465 (2021). 
https://doi.org/10.1038/s41598-021-88802-7

