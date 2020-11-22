# Spatial Prediction via Spatial Weibull Model

![](images/Oil_Production.jpg)

# Contents
1. [Overview](https://github.com/yain22/SWM#1-overview)
2. Research Region: Eagle Ford Shale Reservoir of South Texas
3. Proposed model: Spatial Weibull Model
4. Results of Spatial Prediction at New Test Locations

# 1. Overview

Codes of the paper titled  **"Bayesian hierarchical modeling: application towards production results in the Eagle Ford Shale of South Texas"** are available here. This is a joint project of a Ph.D. student [Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee@stat.tamu.edu) and a Distinguished Professor [Bani K Mallick](https://www.stat.tamu.edu/~bmallick/) (bmallick@stat.tamu.edu) at Texas A&M University. We upload the relevant R codes for the spatial Weibull model for the purpose of the (i) posterior inference (via a Gibbs sampler) to train shale oil wells dataset and (ii) spatial prediction at a new test location based on the trained results. 

The source of data is from [Drillinginfo](https://info.drillinginfo.com/). The data is NOT publically available and there is a cost associated with it charged by the company. Users can buy the data from the website or can get some similar data and will be able to use our codes. Users MUST contact the authors for any use or modification of the uploaded codes for the publication purpose or industrial uses.

**Following R-packages are required:**
1. dplyr
2. tidyr
3. ggplot2
4. maps
5. mvtnorm
6. fields

**Following R-codes are provided:**
1. [SWM.R](https://github.com/yain22/SWM/blob/main/R%20codes/Prediction_SWM.R) : posterior inference; Gibbs sampling algorithm for the  spatial Weibull model
2. [Prediction_SWM.R](https://github.com/yain22/SWM/blob/main/R%20codes/Prediction_SWM.R) : Monte Carlo simulation for the spatial prediction based on the  spatial Weibull model at a new test location
3. [Spatial_Prediction.RMD](https://github.com/yain22/SWM/blob/main/Implementation/Spatial_Prediction.Rmd) : R markdown file to implement the (i) posterior inference and (ii) spatial prediction based on the SWM.R and Prediction_SWM.R.

# 2. Research Region: Eagle Ford Shale Reservoir of South Texas
The Eagle Ford shale reservoir (see ***Figure 1***) is known as possibly the largest single economic development in the history of Texas and ranked as the largest oil and gas development in the world based on capital invested; visit [eaglefordshale.com](https://eaglefordshale.com/) for a detail. In this research, we collected production data from 360 hydraulically fractured horizontal shale oil wells from the Eagle Ford shale reservoir of South Texas; see ***Figure 2*** for the locations of the wells. The time frame of the oil production of 360 well is from January 2011 through June 2017. Hydraulic fracturing horizontal drilling is a drilling technique that makes use of water-based fluids to fracture the reservoir rock formation where the drilling of well takes place in a way that the well runs parallel to the rock formation; See ***Figure 3*** for a schematic example. Completion is defined as the process of making a well ready for the initial production. Completion data (unit) are as follows: (1) first test flowing tubing pressure (psi); (2) ground elevation (ft); (3) measured depth (ft); (4) true vertical depth (ft); (5) upper perforation (ft); (6) lower perforation (ft); (7) perforated interval length (ft); (8) completion count (integer); (9) first test oil volume (bbl); (10) first test oil gravity; and (11) Abstract number. (Unit of abstract is blocks of land units in Texas. It applies to Texas only.)

***Figure 1: Eagle Ford region with three types of petroleum windows. (Source: [United States Energy Information Administration](https://www.eia.gov/))***

![](images/Eagle_Ford_Shale.png)

***Figure 2: Locations of 360 research wells. (Source: [Drillinginfo](https://info.drillinginfo.com/))***

![](images/360_well_locations.JPG)

***Figure 3: A schematic example of a hydraulically fractured horizontal well***

![](images/Hydraulic_Fracturing_explain_detail.png)

# 3. Proposed model: Spatial Weibull Model
The eventual success of petroleum development projects relies on a large degree of well construction costs. As for unconventional reservoirs, because of very low permeability, and a flow mechanism very different from that of conventional reservoirs, estimates for the well construction cost often contain high levels of uncertainty, and oil companies need to make heavy investment in the drilling and completion phase of the wells. We note that the shale oil wells data comprises three aspects of information: oil production rate time series data (denoted by **y**); completion data (denoted by **x**); and well location data (denoted by **s**). ***In this research, our eventual goal is to propose a Bayesian hierarchical model to digest the production results from a shale reservoir region to predict *beforehand* an approximated production behavior of a new well at a new location given specific completion data***. 

We develop a Bayesian hierarchical model, called spatial Weibull model, to exploit the full aspects of the shale oil wells data and fulfill the aforementioned objective. The spatial Weibull model is a fully Bayesian version of non-linear mixed effect model where (i) on the first stage the log-scaled production trajectories from N wells are described by a rate decline curve, Weibull model, and (ii) on the second stage, latent kriging enables the spatial prediction for the production trajectories for a new well at a new location. Additionally, the sparse horseshoe prior on the second stage is used to find some important completion predictors explaining some productivity of wells. See ***Figure 4*** for the hierarchy of the spatial Weibull model and its directed asymmetric graphical model representation.

[SWM.R](https://github.com/yain22/SWM/blob/main/R%20codes/Prediction_SWM.R) is the main R code which implements the Gibbs sampling algorithm for spatial Weibull model to sample from the parameters of the model. Note that the Steps in the code [SWM.R](https://github.com/yain22/SWM/blob/main/R%20codes/Prediction_SWM.R) coincide with the Steps listed in the **Appendix A.2** of the paper. 

***Figure 4: A hierarhcy of the spatial Weibull model (top) and its graphical model representation (bottom)***

![](images/SWM.png)
![](images/graphical_model.png)


# 4. Results of Spatial Prediction at New Test Locations
We use 324 wells as training wells dataset: then, the remaining 36 = 360 - 324 wells are set as the test wells where we can see the performance of the spatial prediction. The locations and production rate trajectories of the training and testing wells are displayed in the ***Figure 5***. 


***Figure 5: Locations (top) and production rate trajectories (bottom) of the training and testing wells***

![](images/Training_Testing_Wells.PNG)
![](images/Production_Traj_Training_Testing_Wells.PNG)


To investigate how the latent kriging of the spatial Weibull model would help improve the predictive performance, we intentionally selected two test wells with API10 numbers set by API10 = 4249332630 (well-1) and API10 = 4231135594 (well-2). See the top panel in ***Figure 6***. Note that the well-1 is positioned near a certain group of training wells, while the well-2 is more isolatedly located, relatively far from other neighboring wells. We assume that if the latent kriging works, then we will have more accurate prediction results for the well-1 than the well-2.

The bottom panel in ***Figure 6*** displays the prediciton results for the two test well. It turns out that the root-mean-square deviations (RMSD) of the two wells based on the spatial Weibull model are 0.178 (well-1) and 0.305 (well-2). (A lower value for the RMSD indicates a better predictive accuracy.) This implies that the production rates of the well-1 has been more accurately predicted than that of the well-2 as we anticipated. From the panels, we can also observe that pointwise posterior predictive 95 % credible interval of the production rates of the well-1 is relatively narrower than that of the well-2 during the production period. To see the effect of the presence of the latent kriging, we also measured RMSDs based on the non-spatial model, which resulted in 0.180 (well-1) and 0.348  (well-2), higher than the RMSDs based on the spatial one, respectively.

***Figure 6: Locations of two test wells (top) and corresponding production results (bottom) with original and log scales***

![](images/Two_Testing_Wells.png)
![](images/Two_Examples.png)

***Table 1*** summarizes prediction results for the estimated ultimate recovery (EUR) for the two test wells: posterior mean, posterior median, and posterior predictive 95 % credible interval of EUR. Because either of the wells did not produce oil for 30 years, there is no reference value to validate a predictive accuracy. Note that the posterior predictive 95 % credible interval of EUR of the well-1 is narrower than that of the well-2.

***Table 1: Summaries of estimated ultimate recoveries for the two test wells***

![](images/Table_EUR.PNG)

Based on the prediction outcomes of the two test wells along with the results from the Subsection 6.1 in the main paper, we conclude that the latent kriging was beneficial in improving the predictive accuracy (the RMSD) and quality (the length of 95 % credible interval). Note that the spatial Weibull model is designed to use the geological proximity information via the geostatistical stochastic processes, hence, this is not a surprising result. More concretely, under the spatial Weibull model, prediction for the well-1 makes use of the proximity advantage, borrowing more information from the nearby wells, while the isolated position of the well-2 renders its prediction sub-optimal. 
