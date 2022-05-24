# temporalECoG
This repository contains Matlab code to reproduce analyses and figures reported in Groen et al., (2022), Temporal dynamics of neural responses in human visual cortex (https://www.biorxiv.org/content/10.1101/2021.08.08.455547v3.full.pdf).

Paper figures can be reproduced by running the scripts located in mkFigures/. 

First, it is necessary to download the BIDS-formatted data from OpenNeuro: 
Visual ECoG dataset https://openneuro.org/datasets/ds004126.

The files located in the derivatives/Groen2022TemporalDynamics/ folder on OpenNeuro should be placed in a directory in this repository named 'analysis' (see function help for analysisRootPath.m). 

The script tde_run.m contains demo code showing how to extract and preprocess the data from the main BIDS directory and how to fit temporal models to the resulting neural time courses. 

The code has dependencies to the following toolboxes:
https://github.com/WinawerLab/ECoG_utils

Some of the data (e.g. retinotopic atlases) were computed using code located in:
https://github.com/BAIRR01/BAIRanalysis

Please contact Iris Groen (i.i.a.groen@uva.nl) for questions.
