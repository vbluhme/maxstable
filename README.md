# Theory and simulation of max-stable processes: With applications to extremal temperatures
Companion repository for master's thesis on max-stable processes written by Viktor Bluhme Jeppesen under the supervision of Thomas Mikosch at the Department of Mathematical Sciences, University of Copenhagen.

**Abstract:**
> In this thesis, we investigate the mathematical properties of max-stable processes and introduce a number of methods for their simulation. We implement simulation by extremal functions and sum-normalized threshold stopping in the statistical programming language `R` using a model-agnostic object-oriented approach. By means of a benchmarking study, we find the extremal functions algorithm to be more efficient than sum-normalization. Using data from the Danish Meteorological Institute, we estimate a max-stable model of Danish summer temperatures. We simulate from this model to estimate return levels of Danish summer temperatures and conditional probabilities of exceedances at high quantiles.

**Overview of files:**
- `sim_maxstable.R`: Object-oriented implementation of extremal functions and sum-normalization algortihms for simulation of max-stable processes, as well as implementations of Brown-Resnick and moving maximum models.
- `fBs_cov.cpp`: C++ helper functions for working with fractional Brownian sheets.
- `figure_spectral_representation.R`: Figure 3.1
- `figure_simulations_sheet.R`: Figures 3.2 and 3.3
- `figure_BR_extcoeff.R`: Figure 3.4
- `figure_naive_BR.R`: Figure 4.1
- `figure_visualize_extremalfuns.R`: Figure 4.2
- `figure_implementation_control.R`: Figures 5.2 and 5.3
- `benchmarking.R`: Benchmarking study; Figures 5.1 and 5.4
- `dmi_1_download.R`: Download and save DMI temperature data
- `dmi_2_estimation.R`: Estimate marginal GEV distribution; normalize temperature data; estimate Brown-Resnick model. Produces all figures and tables in Sections 6.1 and 6.2.
- `dmi_3_simulation.R`: Simulation from max-stable temperature model. Produces all figures in Section 6.3.
