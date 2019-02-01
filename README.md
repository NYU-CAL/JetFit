# JetFit
Gamma-ray Burst Afterglow Light Curve Fitting Tool 

JetFit package can fit the GRB afterglow light curves using boosted fireball model (details see [Wu \& MacFadyen (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aae9de)). It is based on the ScaleFit ([Ryan et al. (2015)](http://iopscience.iop.org/article/10.1088/0004-637X/799/1/3/pdf)). JetFit is in some state between alpha and beta. 

Table is not included in the Github due to the limit of uploading size. Please contact NYU/CAL for a copy and put the table in the proper directory. The current table is the same as the table used in [Wu \& MacFadyen (2018)](https://arxiv.org/abs/1809.06843). We are improving the table by increasing resolution, adding synchrotron absorption and wind circumburst medium. Hopefully, it will come out pretty soon. 

JetFit package consists of three modules: Interpolator, FluxGenerator and Fitter. FluxGenerator can be used separately. 

Two examples are provided to demonstrate the basic use of the package. 
  * Example_FluxGenerator.ipynb shows how to generate synthetic light curves with given parameters. 
  * Example_Fitter.ipynb shows how to fit GW170817 and plot best-fitting light curves and posterior distributions.

