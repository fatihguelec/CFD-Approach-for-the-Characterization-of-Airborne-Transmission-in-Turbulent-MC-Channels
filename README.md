# CFD-Approach-for-the-Characterization-of-Airborne-Transmission-in-Turbulent-MC-Channels

**Overview**

This repository contains MATLAB code and data for modeling airborne pathogen transmission from a molecular communication (MC) perspective using computational fluid dynamics (CFD). The approach simulates the turbulent dispersion of pathogen-laden droplets in the air, which is crucial for understanding how airborne diseases spread, such as COVID-19.

The CFD simulations were carried out in Ansys, and the results were processed and analyzed to statistically characterize airborne pathogen transmission in turbulent environments. This study is based on the statistical characterization of particles received during turbulent flow and the calculation of the infection probability.

**Abstract**

Airborne pathogen transmission mechanisms play a key role in the spread of infectious diseases such as COVID-19. In this work, we propose a computational fluid dynamics (CFD) approach to model and statistically characterize airborne pathogen transmission via pathogen-laden particles in turbulent channels from a molecular communication viewpoint.

To this end, turbulent flows induced by coughing and the turbulent dispersion of droplets and aerosols are modeled using the Reynolds-averaged Navier-Stokes equations, coupled with the realizable k−ε model and the discrete random walk model. Via simulations realized by a CFD simulator, statistical data for the number of received particles are obtained. These data are post-processed to obtain the statistical characterization of the turbulent effect on the reception and to derive the probability of infection.

Our results reveal that turbulence has an irregular effect on the probability of infection, which manifests as a multi-modal distribution—a weighted sum of normal and Weibull distributions. Furthermore, the turbulent MC channel is characterized via multi-modal (i.e., sum of weighted normal distributions) or stable distributions, depending on the air velocity.

**Data Files**

The dataset containing the results from the CFD simulations is available on Zenodo and should be downloaded separately due to its size. After downloading, extract the files and store them in a folder named "Data" in the same directory where the MATLAB code will be run.

Download Dataset: Zenodo DOI: 10.5281/zenodo.13793238
The data consists of .dpm files obtained from CFD simulations in Ansys. These files are essential for reproducing the results and must be placed in the correct directory.

**Code Overview**

The MATLAB code included in this repository corresponds to the figures in the paper. The code performs post-processing of the CFD simulation results to statistically analyze the particle reception and infection probability.

Figures 3, 4, 7-11: The MATLAB code files for each figure can be found in the repository. These files generate the plots and results as presented in the paper.

**How to Use the Code**

Download the dataset from Zenodo and extract it to a folder named "Data" in the same path as the MATLAB code.
Run the MATLAB code corresponding to the figures (Figs. 3, 4, 7-11) to reproduce the results.
Ensure that the necessary Matlab toolboxes for post-processing CFD data are installed.

**Results**

The results demonstrate that turbulent airflow has a significant and irregular effect on the transmission of pathogen-laden particles. The probability of infection in turbulent MC channels can be modeled as a multi-modal distribution or a stable distribution, depending on the specific conditions of the airflow.

Multi-modal distributions are observed at certain air velocities.

These findings provide a deeper understanding of the airborne transmission dynamics in turbulent environments and offer valuable insights into mitigating the spread of airborne infectious diseases.

**Citation**

If you use this repository or the associated dataset in your research, please cite the following paper and dataset:

F. Gulec, F. Dressler and A. W. Eckford, "A Computational Approach for the Characterization of Airborne Pathogen Transmission in Turbulent Molecular Communication Channels," in IEEE Transactions on Molecular, Biological, and Multi-Scale Communications, vol. 9, no. 2, pp. 124-134, June 2023, doi: 10.1109/TMBMC.2023.3273193.

Dataset: Gulec, F., Atakan, B., & Dressler, F. (2022). "CFD Simulation Data for Airborne Pathogen Transmission." Zenodo. DOI: 10.5281/zenodo.13793238
