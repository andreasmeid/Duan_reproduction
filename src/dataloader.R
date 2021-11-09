require(survival); require(pec); require(riskRegression)
require(PWEALL); require(pracma)

#setwd("D:/Duan/hte-prediction-rcts-master/src")
# Load SPRINT or ACCORD datasets.
  load_data <- function(dataset) {
    if (dataset == "sprint") {
      df = read.table("../data/sprint/sprint_cut.csv", header=TRUE, sep=",")
      df["diabetes"] = rep(0, nrow(df))
    } else if (dataset == "accord") {
      df = read.table("../data/accord/accord_cut.csv", header=TRUE, sep=",")
      df["diabetes"] = rep(1, nrow(df))
    } else {
      stop("Unknown dataset")
    }
      y = as.numeric(df[,"cvd"])
      t = as.numeric(df[,"t_cvds"])
      w = as.numeric(df[,"INTENSIVE"])
      df <- df[,-which(names(df) %in% c("cvd", "t_cvds", "INTENSIVE"))]
      return(list(X=df, w=w, y=y, t=t, cols=ncol(df)))
  }
#sprint <- load_data("sprint")
#View(test[["X"]])

  R_searchsorted <- function(a,v) {
    res <- length(a)
    i <- 2
    while(i < length(a)) {
      if (a[i-1] < v & v <= a[i]) {
        res <- i -1
        break
      } else {
        i <- i + 1
      }
      
    }
    return(res)
  }
  
#Calculate inverse propensity of censorship weights for the dataset. 
#Uses a Cox model to model the censoring distribution.
  calculate_ipcw <- function(dataset, time_of_censoring, maxweight=5) {# time_of_censoring zu uebergeben?
    
    df = dataset[["X"]][, 2:(ncol(dataset[["X"]])-1)]
    df$t = floor(as.numeric(dataset[["t"]]))
    df$c = 1 - dataset[["y"]]
  
    cph.fit <- with(df, Surv(t, c))
    #test <- IPCweights(cph.fit)
    coxmodel <- coxph(Surv(t, c)~.,data=df, x=TRUE)
    times <- sort(unique(df$t))
    sf <- predictSurvProb(coxmodel,newdata=dataset[["X"]][, -1],times=times)
    # A matrix with as many rows as NROW(newdata) and as many columns as length(times). 
    # Each entry should be a probability and in rows the values should be decreasing.
    sfmat <- as.matrix(sf);  #dim(sfmat); head(sfmat[1,],50)
    
    weights <- rep(0, nrow(df))
    pb <- txtProgressBar(min=0, max=nrow(df), style=3)
    for (i in 1:nrow(df)) {#
      setTxtProgressBar(pb, value=i)
      sfi <- predictSurvProb(coxmodel,newdata=dataset[["X"]][i,],times=times)
      if (dataset[["y"]][i] ==1) {
        idx <- R_searchsorted(times, dataset[["t"]][i])
        #print(idx)
      } else {
        #print(idx)
        idx <- R_searchsorted(times, time_of_censoring)
      }
      #print(sfmat[i,idx])
      weights[i] <- 1/sfmat[i,idx]
    }
    #[1.24287838 1.19602094 1.1351324  1.23129125 1.31297377 1.23623559

    #weights <- idx <- rep(0, nrow(dataset[["X"]]))
    #for (i in 1:length(weights)) {
    #  if (dataset[["y"]][i] == 1) {
    #    idx[i] = predictSurvProb(coxmodel,newdata=dataset[["X"]][i,],times=dataset[["t"]][i])
    #  } else {
    #    idx[i] = predictSurvProb(coxmodel,newdata=dataset[["X"]][i,],times=floor(time_of_censoring)) # [i]?
    #  }
    #  weights[i] = 1 / idx [i]
    #}
    #weights[weights > maxweight] <- maxweight
    return(weights)
  }

# Take bootstrap sample of the dataset, maintaining proportions as stratified across assignment and outcome.
  bootstrap_dataset <- function(dataset) {
    idxs = 1:nrow(dataset[["X"]])
    treated_cvd = which(dataset[["w"]] == 1 & dataset[["y"]] == 1)
    treated_nocvd = which(dataset[["w"]] == 1 & dataset[["y"]] == 0)
    untreated_cvd = which(dataset[["w"]] == 0 & dataset[["y"]] == 1)
    untreated_nocvd = which(dataset[["w"]] == 0 & dataset[["y"]] == 0)
    idxs = c(
      sample(treated_cvd, size=length(treated_cvd), replace=TRUE),
      sample(treated_nocvd, size=length(treated_nocvd), replace=TRUE),
      sample(untreated_cvd, size=length(untreated_cvd), replace=TRUE),
      sample(untreated_nocvd, size=length(untreated_nocvd),replace=TRUE))
    return(idxs) # oder gesamte Daten dataset[["X"]][idxs,] zurueck?
  }


# Convert a dataset into binary outcomes, with IPCW.
  # (1) Dead < t : dead
  # (2) Alive > t : alive
  # (3) Alive < t : remove from dataset
  # (4) Dead > t : alive
cut_dataset <- function(dataset, cens_time) {
  train = dataset
  idxs = setdiff(1:length(train[["y"]]) ,which(train[["y"]] == 0 & train[["t"]] < cens_time))
  train[["y"]][which(train[["y"]] == 1 & train[["t"]] > cens_time)] = 0
  train[["t"]][which(train[["y"]] == 1 & train[["t"]] > cens_time)] = cens_time
  train_data = list(X=train[["X"]][idxs,],y= train[["y"]][idxs],t= train[["t"]][idxs],w= train[["w"]][idxs])
  ipcw = calculate_ipcw(train_data, time_of_censoring=cens_time)
  train_data = list(X=train[["X"]][idxs,],y= train[["y"]][idxs],t= train[["t"]][idxs],w= train[["w"]][idxs], ipcw=ipcw)
  val_data = list(X=rbind(dataset[["X"]][idxs,],dataset[["X"]][-idxs,]),
                  y= c(dataset[["y"]][idxs],dataset[["y"]][-idxs]),
                  t= c(dataset[["t"]][idxs],dataset[["t"]][-idxs]),
                  w= c(dataset[["w"]][idxs],dataset[["w"]][-idxs]),
                  y_cut= c(train[["y"]][idxs],train[["y"]][-idxs]),
                  cens=c(rep(0, length(idxs)), rep(1, length(dataset[["w"]][-idxs]))))
  return(list(train_data, val_data))
}

combine_datasets <- function(sprint, accord){
  return(list(X=rbind(sprint[["X"]], accord[["X"]]), w=c(sprint[["w"]], accord[["w"]]),
         y=c(sprint[["y"]], accord[["y"]]), t=c(sprint[["t"]], accord[["t"]])))
  }


decision_value_rmst <- function(pred_rr, y, w, t, cens_time, min_km_samples=50) {
  # Return the decision value RMST.
  treat <- which(pred_rr > 0 & w == 1)
  control <- which(pred_rr <= 0 & w == 0)
  if (length(control) < min_km_samples) stop("too few (control) samples") 
  rmsth_1 <-rmsth(y=t[treat],d=y[treat],tcut=floor(cens_time))$rmst
  rmsth_0 <-rmsth(y=t[control],d=y[control],tcut=floor(cens_time))$rmst
  #return(((rmsth_1 * length(w == 1)) + (rmsth_0 * length(w == 0))) / length(y))
  #return(((rmsth_1 * length(treat)) + (rmsth_0 * length(control))) / length(y))
  return(((rmsth_1 * length(treat)) + (rmsth_0 * length(control))) / (length(treat)+length(control)))
}


get_range <- function(scores) {
  lower <- quantile(scores, c(.025)) 
  mean = mean(scores)
  upper = quantile(scores, c(.975)) 
  return(list(lower, mean, upper))
}











### inverse probability of censoring weights
### see van der Laan & Robins (2003)
IPCweights <- function(x, maxweight = 5) {
  
  if (!extends(class(x), "Surv"))
    stop(sQuote("x"), " is not a Surv object")
  
  event <- 1 - x[,2] 
  x[,2] <- 1 - event
  km <- survfit(x ~ 1)
  Ghat <- getsurv(km, times = x[,1]) ## see github issue #54
  Ghat[event == 0] <- 1
  w <- event / Ghat
  w[w > maxweight] <- maxweight
  w
}
### extract survival probabilities
### taken from ipred:::getsurv
### DO NOT TOUCH HERE
getsurv <- function(obj, times)
{
  # get the survival probability for times from KM curve j'
  
  if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
  # <FIXME: methods may have problems with that>
  class(obj) <- NULL
  # </FIXME>
  lt <- length(times)
  nsurv <- times
  
  # if the times are the same, return the km-curve
  
  if(length(times) == length(obj$time)) {
    if (all(times == obj$time)) return(obj$surv)
  }
  
  # otherwise get the km-value for every element of times separatly
  
  inside <- times %in% obj$time
  for (i in (1:lt)) {
    if (inside[i])
      nsurv[i] <- obj$surv[obj$time == times[i]]
    else  {
      less <- obj$time[obj$time < times[i]]
      if (length(less) == 0)
        nsurv[i] <- 1
      else
        nsurv[i] <- obj$surv[obj$time == max(less)]
    }
  }
  nsurv
}
