README

Code accompanying:

Fang and Willett (2026)
Climatic cycles and potential bias in the relationship between channel steepness and cosmogenic erosion rates

Version 1.0

==================================================
Overview
==================================================

This repository contains the Python code and supporting data used for the numerical experiments and natural case-study analyses presented in the accompanying manuscript.

The repository includes two main components:

a one-dimensional river-profile model used to investigate how cyclic precipitation affects erosion-rate and channel-steepness metrics;
a set of scripts used to calculate the Cosmogenic-constrained cycle-mean erosion rate, Ecos-Qm, and reproduce the analysis of granite catchments in the Three Rivers region.

The one-dimensional model simulates sinusoidal precipitation cycles with different combinations of precipitation amplitude(Pa) and mean precipitation (Pm).

The case study scripts combine catchment data with palaeoprecipitation information to calculate Ecos-Qm.

==================================================
Repository contents
==================================================

One-dimensional river-profile model

1.climatic_cycles_1d.py
	Runs the one-dimensional longitudinal river-profile model used in the manuscript. The script simulates river-profile evolution under sinusoidal precipitation cycles; allows different precipitation amplitudes and mean precipitation rates to be tested; calculates Ein, Ecos, Ecos-Qm, ks, ks-Q, ks-Qeff


Three Rivers region case study

2.climatic_cycle_preparation.py
	Prepares the precipitation and catchment data required for the Three Rivers region case study. It reads the relevant input datasets, processes the palaeoprecipitation records, and generates the intermediate input files used by the main analysis script.

3.climatic_cycle_functions.py
	Contains the functions used by the case study scripts. This file is imported by the other scripts and does not need to be executed separately.

4.climatic_cycle_main.py
	Runs the main Three Rivers region case study analysis. The script calculates Ecos-Qm and associated channel-steepness metrics for the study catchments and generates the case study results and figures presented in the manuscript.

==================================================
Supporting data
==================================================

The following supporting files are included in this package.

granite_basin.shp
	Polygon shapefile of the study catchments.

granite_basin.xlsx
	Dataset containing production rates, cosmogenic erosion rates, ksn-Q, and ksn-Qeff, calculated following the methods described in the accompanying manuscript.

proxy_benthic_Lisiecki2005.xls
	Benthic oxygen-isotope proxy record derived from Lisiecki and Raymo (2005).

proxy_insolation_Laskar2004_10Myr.xls
	Orbital-insolation record derived from Laskar et al. (2004), covering the past 10 million years.

==================================================
External data
==================================================

The CHELSA-TraCE21k precipitation dataset is not included because of its large file size.

It can be downloaded from:
https://www.chelsa-climate.org/models/chelsa-trace21k


==================================================
Citation
==================================================

If you use this code, please cite:

Fang and Willett (2026).
<Climatic cycles and potential bias in the relationship between channel steepness and cosmogenic erosion rates>. (submitted to Esurf)
Contact Xianjun Fang (xianjun.fang@eaps.ethz.ch)

==================================================
References
==================================================

Karger, D.N., Nobis, M.P., Normand, S., Graham, C.H., Zimmermann, N.E., 2021. CHELSA-TraCE21k v1.0. Downscaled transient temperature and precipitation data since the last glacial maximum. https://doi.org/10.5194/cp-2021-30
Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B., 2004. A long-term numerical solution for the insolation quantities of the Earth. Astron. Astrophys. 428, 261–285. https://doi.org/10.1051/0004-6361:20041335
Lisiecki, L.E., Raymo, M.E., 2005. A Pliocene‐Pleistocene stack of 57 globally distributed benthic δ 18 O records. Paleoceanography 20, 2004PA001071. https://doi.org/10.1029/2004PA001071


The Python source code is licensed under the MIT License. The original catchment data produced for this study are licensed under the Creative Commons Attribution 4.0 International License. Third-party data remain subject to the terms of their original providers.
