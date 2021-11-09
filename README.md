### Reproducing Duan et al. "Clinical Value of Predicting Individual Treatment
Effects for Intensive Blood Pressure Therapy" Circ Cardiovasc Qual Outcomes 2019;12:e005010.

Last update: July 2019.

---

Code after translating the original analysis code in Python to R from the paper:

[Clinical value of predicting individual treatment effects for intensive blood pressure therapy](https://www.ahajournals.org/doi/10.1161/CIRCOUTCOMES.118.005010).<br/>
Tony Duan, Pranav Rajpurkar, Dillon Laird, Andrew Y. Ng, Sanjay Basu.
_Circulation: CQO_, 2019

#### Model development

The predictive modeling is based on the X-learner algorithm stored in the lib-subfolder.

#### Evaluation

According to the original publication, evaluation relies on R implementations of:

1. C-statistic-for-benefit 
2. Decision value of restricted mean survival time (RMST) 
3. Comparisons of absolute risk reductions in buckets of patients recommended an intensified blood pressure lowering or not.

For all statistics, bootstrap confidence intervals were calculated by resampling of the dataset.

#### Reproduction

The data required to run the code relates to the SPRINT and ACCORD-BP datasets to be downloaded from [BioLINCC] (https://biolincc.nhlbi.nih.gov/home/).

The make.R file runs the analyses with the default 250 bootstrap samples. 

