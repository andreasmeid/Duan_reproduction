# Framingham risk functions

framingham=function(Sex=c("M","F"), Age, BMI, Smoke=c("Y","N"), DM=c("Y","N"),SBP, Meds=c("Y","N"), Totalchol, HDL , method=c("Office","Full"))
{
  #  John W Pickering 2016, October 2016. john.pickering @ otago.ac.nz (c) Creative Commons Attribution-Share Alike 4.0 International
  #  use at own risk (I checked against some data I had, but not guaranteed)
  #  See: https://www.framinghamheartstudy.org/risk-functions/cardiovascular-disease/10-year-risk.php
  #  Note:  cholesterol is in mg/dL.  To convert from mmol/L divide by 0.02586
  
  # Converstions
  age <- log(Age)
  bmi <- log(BMI)
  sbp <- log(SBP)  
  
  sex = ifelse(Sex=="M",1,0)
  smoke <- ifelse(Smoke=="Y",1,0)
  dm <- ifelse(DM=="Y",1,0)
  meds <- ifelse(Meds=="Y",1,0)
  nomeds <- ifelse(Meds=="N",1,0)
  
  sbp_meds <- sbp*meds
  sbp_nomeds<-sbp*nomeds
  
  df <- data.frame(sex,age,bmi,smoke,sbp_meds,sbp_nomeds,dm)
  n <- dim(df)[1]
  r<-matrix(nrow=n,ncol=1)
  
  # Method and calculation
  
  if(method=="Full")
  {
    df$totalchol <- log(Totalchol)
    df$hdl <- log(HDL)
    attach(df)
    for (i in 1: n)
    {
      ifelse(sex[i]==0,
             r[i] <- (1-0.95012^(exp(((2.32888*age[i])+(1.20904*totalchol[i])+(-0.70833*hdl[i])+(2.76157*sbp_nomeds[i])+(2.82263*sbp_meds[i])+(0.52873*smoke[i])+(0.69154*dm[i]))-26.1931)))*100,
             r[i] <- (1-0.88936^(exp(((3.06117*age[i])+(1.12370*totalchol[i])+(-0.93263*hdl[i])+(1.93303*sbp_nomeds[i])+(1.99881*sbp_meds[i])+(.65451*smoke[i])+(.57367*dm[i]))-23.9802)))*100 )
    }
  }  
  
  if(method=="Office")
  {
    attach(df)
    for (i in 1: n)
    {
      ifelse(sex[i]==0,
             r[i] <- (1-0.94833^(exp(((2.72107*age[i])+(0.51125*bmi[i])+(2.81291*sbp_nomeds[i])+(2.88267*sbp_meds[i])+(0.61868*smoke[i])+(.77763*dm[i]))-26.0145)))*100,
             r[i] <- (1-0.88431^(exp(((3.11296*age[i])+(.79277*bmi[i])+(1.85508*sbp_nomeds[i])+(1.92672*sbp_meds[i])+(.70953*smoke[i])+(.53160*dm[i]))-23.9388)))*100  )
    }
  }  
  return(r)
}



BETA_WOMEN <- c(2.72107, # log age
                0.51125,          # log BMI
                2.81291,          # log SBP (not treated)
                2.88267,          # log SBP (treated)
                0.61868,          # smoking
                0.77763)          # diabetes

BETA_MEN <- c(3.11296,0.79277,1.85508,1.92672,0.70953,0.53160)

S_WOMEN <- c(0.99781, 0.99107, 0.98756, 0.98145, 0.97834, 0.97299, 0.96635, 0.96210, 0.95580, 0.94833)
S_MEN <- c(0.99432, 0.98281, 0.97433, 0.95748, 0.94803, 0.93434, 0.92360, 0.91482, 0.90342, 0.88431)

CONST_WOMEN <- 26.0145
CONST_MEN <- 23.9388

calc_frs <- function(X, b, surv, const) {
  #Simple Non-Laboratory Framingham Risk Score (FRS) Calculation.
  
  #Parameters
  #----------
  #  X : array or list
  #Input variables for log-age, log-bmi, log sbp (not treated),
  #log sbp (treated), smoking, diabetes
  #b : array or list
  #Variable coefficients
  #surv : float  #Baseline survival
  #const : float
  #Model intercept
  return(1 - surv^(exp(X%*%(b) - const))) # ** power exponent
}

frs <- function(gender, time, age, bmi, sbp, ht_treat, smk, dia) {
  #10-year risk for women, calculated using the Simple Non-Laboratory Framingham Risk Score (FRS) Calculation.
  #Parameters
  #----------
  #gender : char # Gender of subject ('M' or 'F')
  #time : int #Time horizon for risk calculation. Must be between 1 and 10.
  #age : numeric #Age of subject
  #bmi : numeric # BMI of subject
  #sbp : numeric #Systolic blood pressure of subject
  #ht_treat : bool or int #Treatment for hypertension (True or False)
  #smk : bool or int # Subject is smoker (True or False)
  #dia : bool or int # Subject has diabetes (True or False)
  if (time<1 | time>10) stop('Risk can only be calculated for 1 to 10 year time horizon')
  ht = as.logical(ht_treat)
  X = c(log(age), log(bmi), log(sbp)*(1-ht), log(sbp)*ht, as.logical(smk), as.logical(dia))
  if (gender=='F') {
    ret <- calc_frs(X, BETA_WOMEN, S_WOMEN[time-1], CONST_WOMEN)
  } else if (gender=='M') {
    ret <- calc_frs(X, BETA_MEN, S_MEN[time-1], CONST_MEN)
  } else {
    stop('Gender must be specified as M or F')
  }
  #test
    #time=10
    #X = c(log(35), log(24.3), log(122)*(1-FALSE), log(122)*FALSE, TRUE, FALSE)
    #calc_frs(X, BETA_WOMEN, S_WOMEN[time-1], CONST_WOMEN)
  #0.029352227213368165
  #frs(gender='F', time=10, age=35, bmi=24.3, sbp=122, ht_treat=False, smk=True, dia=False)
  return(ret)
}


