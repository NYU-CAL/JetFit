# JetFit
Gamma-ray Burst Afterglow Light Curve Fitting Tool. 

JetFit package can fit the GRB afterglow light curves using boosted fireball model (details see [Wu \& MacFadyen (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aae9de)). It is based on the ScaleFit ([Ryan et al. (2015)](http://iopscience.iop.org/article/10.1088/0004-637X/799/1/3/pdf)). JetFit is in some state between alpha and beta. 

`Table.h5` contains the characteristic spectral functions, which are used to generate synthetic light curves. 
The table is almost the same as the table used in [Wu \& MacFadyen (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aae9de). We are improving the table by increasing resolution, adding synchrotron absorption and wind circumburst medium. Hopefully, it will come out pretty soon. 

# Installation
Anaconda is recommened to be installed. The detailed python package list can be found in package_list.txt. To set up the exactly the same virtual environment, run ```conda create --name JetFit --file package_list.txt``` 


to create a environment `JetFit`. Then, run `conda activate JetFit` to activate the environment. It is all set. 

JetFit package consists of three classes: Interpolator, FluxGenerator and Fitter. FluxGenerator can be used separately. 

Two examples are provided to demonstrate the basic use of the package. 
  * Example_FluxGenerator.ipynb shows how to generate synthetic light curves with given parameters. 
  * Example_Fitter.ipynb shows how to fit GW170817 and plot best-fitting light curves and posterior distributions.
