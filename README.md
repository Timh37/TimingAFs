

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
* GESLA3 data (download [here](https://gesla787883612.wordpress.com/downloads/))
* FLOPROS data from Tiggeloven et al. (2020) (download [here](https://nhess.copernicus.org/articles/20/1025/2020/))
* AR6 Sea-level projections at GESLA3 locations (uploaded [here]())

<p align="right">(<a href="#top">back to top</a>)</p>

## Workflow
1. Derive daily maxima from GESLA3 data and write them to CSV files using the script [daily_max_GESLA3_to_csv.py](https://github.com/Timh37/TimingAFs/blob/main/GPD_analysis/daily_max_GESLA3_to_csv.py).
2. Load the CSV files and apply the method of Solari et al. (2017) to the daily maxima using [gpdfit_solaris_thres_gesla3.R](https://github.com/Timh37/TimingAFs/blob/main/GPD_analysis/gpdfit_solaris_thres_gesla3.R). This script automatically selects the extreme threshold and outputs central-estimate general Pareto distribution parameters and their uncertainty to NetCDF files.
3. To project the timing of AFs, [project_AF_timing_ar6wfs.py](https://github.com/Timh37/TimingAFs/blob/main/project_timing/project_AF_timing_ar6wfs.py)  takes the GPD parameter estimates and fetches the FLOPROS estimates nearest to the GESLA3 locations. Using the GPD fits, the return curves and required SLR is computed. The timing of AFs is then projected by interpolating the timing of projected SLR onto required SLR, and all is stored in a single NetCDF file for all locations.

<p align="right">(<a href="#top">back to top</a>)</p>

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
We thank Sebastian Solari for sharing his automatic threshold selection code and Timothy Tiggeloven for elucidating the FLOPROS estimates. R.E.K. and M.O. were supported by the National Science Foundation (NSF) as part of the Megalopolitan Coastal Transformation Hub (MACH) under NSF award ICER-2103754. V.M.-S. and A.B.A.S. were supported by PROTECT. This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 869304, PROTECT contribution number XX.

<p align="right">(<a href="#top">back to top</a>)</p>

## Contact
Tim Hermans
t.h.j.hermans@uu.nl

<p align="right">(<a href="#top">back to top</a>)</p>

