# DMD & HDMD of the Annular Modes
This repository contains scripts/codes to calculate Dynamic mode decomposition (DMD) and Hankel Dynamic Mode Decomposition (HDMD) of the Annular Modes based on Matlab and NCL

1. dmd_u_erainterim.m includes the codes to calculate dynamic mode decomposition of the zonal mean zonal wind anomalies (p, lat) from ERA-Interim.

2. eof_u_erainterim.ncl includes the codes to calculate empirical orthogonal functions of the zonal mean zonal wind anomalies from ERA-Interim.

3. corr_dmd_u_erainterim.ncl includes the codes to find pairs of DMD modes that are highly correlated with the EOF modes.

4. dmd_u_erainterim_plot.ncl includes the codes to plot the best pairs of DMD modes.

<p align="center">
  <img src="https://github.com/sandrolubis/DMD_HDMD_Annular_Modes/blob/main/example/dmd_ref_new_crop.png" width="800">
</p>


**Figure 1.** Principal oscillation patterns (POPs) or dynamical mode decomposition (DMD) patterns (color shading) superimposed with the two leading EOF modes (contour lines) for the Southern Hemisphere zonal mean large-scale circulation. (a) Real and (b) imaginary parts of the dominant Southern Hemisphere POPs from year-round ERA-Interim data, calculated at θ =11 days (representing the modes onto which the first two EOFs project the most strongly). The corresponding eigenvalue is λ = -0.0974 ± 0.0425i day−1. The contour lines denote the EOF1 (EOF2) pattern superimposed with the real (imaginary) parts of the POPs (**Lubis and Hassanzadeh, in prep**).


5. hdmd_u_erainterim.m includes the codes to calculate Hankel dynamic mode decomposition of the zonal mean zonal wind anomalies (p, lat) from ERA-Interim.

6. eof_u_erainterim.ncl includes the codes to calculate empirical orthogonal functions of the zonal mean zonal wind anomalies from ERA-Interim.

7. corr_hdmd_u_erainterim.ncl includes the codes to find pairs of HDMD modes that are highly correlated with the EOF modes.

8. hdmd_u_erainterim_plot.ncl includes the codes to plot the best pairs of HDMD modes.

9. GCM_HDMD.m is alculating the basis functions (EOFs or DMD modes or Hankel-DMD modes depending on the entries)
