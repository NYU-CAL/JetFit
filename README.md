# JetFit
Gamma-ray Burst Afterglow Light Curve Fitting Tool.

JetFit package can fit the GRB afterglow light curves using boosted fireball model (details see [Wu \& MacFadyen (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aae9de)). It is based on the ScaleFit ([Ryan et al. (2015)](http://iopscience.iop.org/article/10.1088/0004-637X/799/1/3/pdf)). JetFit is in some state between alpha and beta.

`Table.h5` contains the characteristic spectral functions, which are used to generate synthetic light curves.
The table is almost the same as the table used in [Wu \& MacFadyen (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aae9de). We are improving the table by increasing resolution, adding synchrotron absorption and wind circumburst medium. Hopefully, it will come out pretty soon.

# Installation
Anaconda is recommened to be installed. The detailed python package list can be found in package_list.txt. To set up the virtual environment, run command `conda create --name JetFit --file package_list.txt` in terminal. Then, run `conda activate JetFit` to activate the environment.

# Usage
To fit the light curve from GW170817, run `python Example_Fitter.py`.

# Module Description
JetFit package consists of three classes: Interpolator, FluxGenerator and Fitter. FluxGenerator can be used separately.

 ## InterpolatorClass:
  * `_LoadTable`: load characteristic spectral function table.
  * `_SetScale`: set scales for the table. By default, f_peak and tau are in log scale.
  * `_GetInterpolator`: use `scipy.interpolate.RegularGridInterpolator` as interpolator.
  * `GetValue`: get values for characteristic spectral function at specific position (tau, Eta0, GammaB, theta_obs).

 ## FluxGeneratorClass:
  * `GetTaus`: rescale the observational time in second.
  * `GetTransformedValue`: get the transformed values F_peak, nu_c and nu_m (Ryan et al. 2015).
  * `GetSpectral`: calculate synthetic light curves F_nu (Sari et al. 1998).
  * `GetIntegratedFlux`: calculate synthetic integrated light curves.

 ## ScaleFitClass:
  * `_SetFitParameter`: set up ScaleFit, e.g. parameter scales and parameter bounds.
  * `LoadData`: load observational data.
  * `GetSampler`: initialize sampler from  emcee package.
  * `BurnIn`: burn in and prepare for sampling.
  * `RunSampler`: run sampling procedure and save results to local drive.
