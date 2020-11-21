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
**SWM.R** is the main file which implements the Gibbs sampler for the Spatial Weibull Model (SWM) used to learn parameters from the training wells. Steps in the Gibbs sampler are coincided with the steps listed in the **Appendix A.2** of the paper.

