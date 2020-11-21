# Spatial Prediction via Spatial Weibull Model

Codes of the paper titled  **"Bayesian hierarchical modeling: application towards production results in the Eagle Ford Shale of South Texas"** are avaiable here. This is a joint project of a Ph.D. student [Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee@stat.tamu.edu) and a Distinguished Professor [Bani K Mallick](https://www.stat.tamu.edu/~bmallick/) (bmallick@stat.tamu.edu) at Texas A&M University. We upload the relevant R codes for the spatial Weibull model (SWM) for the purpose of the (i) posterior inference (via a Gibbs sampler) to learn from training data and (ii) spatial prediction at a new test location. 

Source of data is from [Drillinginfo](https://info.drillinginfo.com/). The data is NOT publically available and there is a cost associated to it charged by the company. Users can buy the data from the website or can get some similar data and will be able to use our codes. Users MUST contact the authors for any use or modification for our codes for the publication purpose or industrial uses.

**Following R-packages are required:**
1. dplyr
2. tidyr
3. ggplot2
4. maps
5. mvtnorm
6. fields

**Following R-codes are provided:**
1. SWM.R : posterior inference; Gibbs sampling algorithm for the SWM
2. Prediction_SWM.R : Monte Carlo simulation for the spatial prediciton based on the SWM at a new test location
3. Spatial_Prediction.RMD : R markdown file to implement the (i) posterior inference and (ii) spatial prediction based on the SWM.R and Prediction_SWM.R.

# 1. Research Region: Eagle Ford Shale Reservoir of South Texas
We research 360 hydraulically-fractured horizontal shale oil wells data collected from the Eagle Ford shale reservoir of South Texas. The Eagle Ford shale reservoir (see *Figure 1* is known as possibly the largest single economic development in the history of Texas and ranked as the largest oil and gas development in the world based on capital invested; visit [eaglefordshale.com](https://eaglefordshale.com/) for a detail. 

Time frame of the oil production of 360 well is from January 2011 through June 2017. Hydraulic-fracturing-horizontal drilling is a well completion technique which makes use of water-based fluids to fracture the reservoir rock formation where the drilling of well takes place in a way that the well runs parallel to the rock formation; See *Figure 2* for a schematic example of a hydraulically fractured horizontal well. Completion is defined as the process of making a well ready for the initial production.  

*Figure 1: Eagle Ford region with three types of petroleum windows. (Source: United States Energy Information Administration)*

![](images/Eagle_Ford_Shale.png)

*Figure 2: A schematic example of a hydraulically fractured horizontal well*

![](images/Hydraulic_Fracturing_explain_detail.png)



## Spatial Weibull Model
Spatial Weibull Model (SWM) 

**SWM.R** is the main file which implements the Gibbs sampling algorithm for the Spatial Weibull Model (SWM) to sample from the parameters the model. Note that the Steps in the code **SWM.R** are coincided with the Steps listed in the **Appendix A.2** of the paper. 

![](images/graphical_model.png)
