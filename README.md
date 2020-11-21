# Spatial Prediction via Spatial Weibull Model

Codes of the paper titled  **"Bayesian hierarchical modeling: application towards production results in the Eagle Ford Shale of South Texas"** are avaiable here. This is a joint project of [Se Yoon Lee](https://sites.google.com/view/seyoonlee) and [Bani K Mallick](https://www.stat.tamu.edu/~bmallick/). We upload the relevant R codes for the spatial Weibull model (SWM) for the purpose of the (i) posterior inference and (ii) spatial prediction at a new test location. 

**Following R-packages are required:**
1. dplyr
2. tidyr
3. ggplot2
4. maps
5. mvtnorm
6. fields

## Research Region: Eagle Ford Shale Reservoir of South Texas
We research 360 hydraulically-fractured horizontal oil wells data collected from the Eagle Ford shale reservoir of South Texas. Time frame of the oil production of 360 well is from January 2011 through June 2017. Hydraulic-fracturing-horizontal drilling is a well completion technique which makes use of water-based fluids to fracture the reservoir rock formation where the drilling of well takes place in a way that the well runs parallel to the rock formation. Completion is defined as the process of making a well ready for the initial production, and completion data is defined as any meaningful data involved in completion procedure. Units of completion data are psi (pounds per square inch), ft (feet), bbl (barrel), etc. Unit of location of well is (longitude, latitude) such that decimal degrees has been used in WGS84 coordinate reference system. Unit of oil production time series data is bbl/month (barrel per month). 

![](images/Eagle_Ford_Shale.png)
![](images/Hydraulic_Fracturing_explain_detail.png)



## Spatial Weibull Model
Spatial Weibull Model (SWM) 

**SWM.R** is the main file which implements the Gibbs sampling algorithm for the Spatial Weibull Model (SWM) to sample from the parameters the model. Note that the Steps in the code **SWM.R** are coincided with the Steps listed in the **Appendix A.2** of the paper. 

![](images/graphical_model.png)
