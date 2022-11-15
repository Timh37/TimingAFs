

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#TimingAFs">TimingAFs</a></li>
    <li><a href="#Prerequisites-and-built-with">Prerequisites and built with</a></li>
    <li><a href="#Acknowledgements">Acknowledgements</a></li>
    <li><a href="#Contact">Contact</a></li>
    
  </ol>
</details>

# TimingAFs
TimingAFs is a novel framework for projecting the timing of frequency amplifications of extreme sea levels due to sea-level rise. The framework allows the user to compute when a certain amplification of a benchmark frequency will occur, instead of how much it will amplify by 2100 or another arbitrary future year. Instead of using the historical centennial event as a benchmark for amplification factors, estimated flood protection standard from Tiggeloven et al. (2020) are used as a benchmark frequency. The projections therefore indicate the timing of decreases in the degree of coastal flood protection at a given location. For the projections, the sea-level projections of IPCC AR6 are used. The extreme sea level distributions are derived from GESLA3 tide gauge data using a peak-over-threshold analysis with an automatic threshold selection based on Solari et al. (2017).

The code in this repository was used for the manuscript:
```
The Timing of Decreasing Coastal Protection Levels Due to Sea-Level Rise (under review), Hermans et al.
```
Please cite it when using this repository for your own studies.

<p align="right">(<a href="#top">back to top</a>)</p>

## Required input data
* GESLA3 data ([download here](https://gesla787883612.wordpress.com/downloads/))
* FLOPROS data from [Tiggeloven et al. (2020)](https://nhess.copernicus.org/articles/20/1025/2020/)
* AR6 Sea-level projections (full sample total relative projections at GESLA3 data can be found here: to-do)

## Workflow
1. Derive daily maxima from GESLA3 data using the scripts in */TimingAFs/GPD_analysis/*.
2. Apply the method of Solari et al. (2017) to the daily maxima to automatically select the extremes threshold at each GESLA3 site, and derive best-estimate general Pareto distribution parameters and their uncertainty using bootstrapping with the script */TimingAFs/GPD_analysis/gpdfit_solaris_thres_gesla3.R*. This script reads in daily maxima stored to csv files for each site.
3. Use the resulting GPD distributions to generate return curves, fetch benchmark frequencies from FLOPROS, compute required SLR and interpolate the timing of projected SLR to required SLR to project the timing of AFs, using the scripts in */TimingAFs/project_timing/*

## Prerequisites and built with

* Python
  * numpy
  * matplotlib
  * xarray
  * scipy
  * os
  * cartopy
  * cmocean
  * shapely
  * gesla

<p align="right">(<a href="#top">back to top</a>)</p>

## Acknowledgements

<p align="right">(<a href="#top">back to top</a>)</p>

## Contact
Tim Hermans
t.h.j.hermans@uu.nl

<p align="right">(<a href="#top">back to top</a>)</p>

