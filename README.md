# Spatial Prediction via Spatial Weibull Model

Codes of the paper titled  **"Bayesian hierarchical modeling: application towards production results in the Eagle Ford Shale of South Texas"** are avaiable here. This is a joint project of **Se Yoon Lee** and  **Bani K Mallick**. 

**Following R-packages are required:**
1. dplyr
2. tidyr
3. ggplot2
4. maps
5. mvtnorm
6. fields

## Gibbs Sampler of the Spatial Weibull Model
**SWM.R** is the main file which implements the Gibbs sampling algorithm for the Spatial Weibull Model (SWM) to sample from the parameters the model. Note that the Steps in the code **SWM.R** are coincided with the Steps listed in the **Appendix A.2** of the paper. 

![](images/graphical_model.png)
