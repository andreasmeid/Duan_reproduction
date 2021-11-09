require(PWEALL); require(pracma); 
require(Matching); require(Hmisc);
require(binr)
library(data.table)

# Return the empirical 95% range of a statistic.
  get_range <- function(scores) {
    lower <- quantile(scores, c(.025)) 
    mean = mean(scores)
    upper = quantile(scores, c(.975)) 
    return(list(lower, mean, upper))
  }
  # test <- rnorm(20,1,1)
  

# Return p-values for the survival rate in each of the buckets.
  wald_test <- function (pred_rr=pred_rr[cens_ind], y=y_cut[cens_ind], w=w[cens_ind]) {
  # buckets of benefit (absolute risk reduction [ARR] >0) and no benefit (ARR ???0)
    selection_criterion1 <- which(as.numeric(pred_rr) > 0)
    selection_criterion2 <- which(as.numeric(pred_rr) <= 0)
    bucket_ben <- data.frame(#pred_rr=pred_rr[selection_criterion1],
                             y=y[selection_criterion1], w=w[selection_criterion1])  
    bucket_noben <- data.frame(#pred_rr=pred_rr[selection_criterion2],
                             y=y[selection_criterion2], w=w[selection_criterion2])  
    bucket_ben$w <- ifelse(bucket_ben$w == 1, "y", "n"); #table(bucket_ben$y, bucket_ben$w)
    p_ben <- chisq.test(bucket_ben$w, bucket_ben$y, correct=FALSE)[["p.value"]]
    
    bucket_noben$w <- ifelse(bucket_noben$w == 1, "y", "n"); #table(bucket_noben$y, bucket_noben$w)
    p_noben <- chisq.test(bucket_noben$w, bucket_noben$y, correct=FALSE)[["p.value"]]

    return(list(p_ben, p_noben))
  }

# Discrimination was further assessed by bucketing participants into subgroups of predicted benefit 
  # (those with predicted ARR suggesting lower CVD risk with intensive treatment) and no benefit 
  # (those with zero or negative predicted ARR).
  
bucket_arr <- function(pred_rr, y, w) {
  # Evaluation of model bucketing.
  tempdf <- data.frame(pred_rr=pred_rr, y=y, w=w)
  
  
  arr_ben = mean(tempdf[which(tempdf$pred_rr > 0 & tempdf$w == 0),]$y) -  mean(tempdf[which(tempdf$pred_rr > 0 & tempdf$w == 1),]$y)
  arr_noben = mean(tempdf[which(tempdf$pred_rr <= 0 & tempdf$w == 0),]$y) -  mean(tempdf[which(tempdf$pred_rr <= 0 & tempdf$w == 1),]$y)
  return(list(arr_ben, arr_noben))
}

# decision_value_rmst moved to dataloader.R


c_statistic <- function(c_pred_rr, c_y, c_w) {
  #c_pred_rr <- cens_pred_rr[idxs]
  #c_y <- cens_y_cut[idxs]
  #c_w <- cens_w[idxs]
  # Return concordance-for-benefit, the proportion of all matched pairs with
  # unequal observed benefit, in which the patient pair receiving greater
  # treatment benefit was predicted to do so.
  
  # ensure results are reproducible
    #set.seed(123)
    #assert len(pred_rr) == len(w) == len(y)
  
  # match all pairs on predicted benefit
    tuples = data.frame(c_pred_rr=c_pred_rr, c_y=c_y, c_w=c_w)
    untreated = tuples[which(tuples$c_w == 0),]
    treated = tuples[which(tuples$c_w == 1),]
  
  # randomly subsample to ensure every person is matched
    if (nrow(treated) < nrow(untreated)) {
      untreated <- untreated[sample(1:nrow(untreated), size=nrow(treated), replace=FALSE),]
    } else if (nrow(untreated) < nrow(treated)) {
      treated <- treated[sample(1:nrow(treated), size=nrow(untreated), replace=FALSE),]
    }
    untreated <- untreated[order(untreated$c_pred_rr),] #untreated = sorted(untreated, key=lambda t: t[0])
    treated <- treated[order(treated$c_pred_rr),] # treated = sorted(treated, key=lambda t: t[0])
    
    make_pairs <- function(var1, var2) {
      return(paste(var1, var2, sep="-"))
    }
    
    obs_benefit_dict <- function(paar) {
      res <- rep(0,length(paar))
      for (i in 1:length(paar)) {
        if (paar[i] == "0-0" | paar[i] == "1-1") {
          res[i] <- 0
        } else if (paar[i] == "0-1") {
          res[i] <- -1
        } else if (paar[i] == "1-0") {
          res[i] <- 1
        }
      }
      
      return(res)
    }
    
    obs_pairs <- make_pairs(untreated$c_y, treated$c_y)
    
  # calculate observed and predicted benefit for each pair
    
    obs_benefit = obs_benefit_dict(obs_pairs)
    pred_benefit = rowMeans(cbind(untreated$c_pred_rr, treated$c_pred_rr))
  
    cindex <- rcorr.cens(pred_benefit, obs_benefit)
    c.benefit <- cindex["C Index"][[1]]
    
  # iterate through all (N choose 2) pairs
    #count <- total <- 0
    #for (i in 1:(length(obs_pairs)-1)) {
    #  for (j in (i + 1):length(obs_pairs)) {
    #    if (obs_benefit[i] != obs_benefit[j]) {
    #      if (((obs_benefit[i] < obs_benefit[j]) &
    #           (pred_benefit[i] < pred_benefit[j])) |
    #          ((obs_benefit[i] > obs_benefit[j]) &
    #           (pred_benefit[i] > pred_benefit[j]))) {
    #        count <- count + 1
    #      } 
    #    }
    # }
    #}
  #return(count / length(obs_pairs))
  return(c.benefit)
}

c_for_benefit <- function(pred.ben=static.prob, treatment=val_dat2$DOAK, y=val_dat2$y, 
                          T_A="Dabi", T_B="Apixa", covs=NULL,#covs <- c("val_dat2$HighAge", "val_dat2$Female", "val_dat2$TE" )
                          matching_for="benefit", MatchingPlot= TRUE) {
  
  if (matching_for != "benefit" & is.null(covs)==FALSE) { # Matching on covariates
    eval(parse(text=paste("covdata <- data.frame(pred.ben, ",  paste(covs, collapse=','), ")", sep="")))
    set.seed(1) # For reproducibility
    rr <- Match(Tr=treatment==T_B, X=covdata, M=1,ties=F,replace=FALSE)
    ind.A <- rr$index.control
    ind.B <- rr$index.treated
    ### Calculation of predicted and observed benefit in matched pairs
    pred.ben.A <- pred.ben[ind.A]
    pred.ben.B <- pred.ben[ind.B]
    pred.ben.avg <- (pred.ben.A+pred.ben.B)/2
    obs.out.A <- y[ind.A]
    obs.out.B <- y[ind.B]
    obs.ben <- obs.out.A-obs.out.B
    
  } else { # matching on benefit
    ind.A <- which(treatment==T_A)
    order.A <- order(pred.ben[ind.A])
    ind.A <- ind.A[order.A]
    ind.B <- which(treatment==T_B)
    order.B <- order(pred.ben[ind.B])
    ind.B <- ind.B[order.B]
    ### Calculation of predicted and observed benefit in matched pairs
    pred.ben.A <- pred.ben[ind.A]
    pred.ben.B <- pred.ben[ind.B]
    pred.ben.avg <- (pred.ben.A+pred.ben.B)/2
    obs.out.A <- y[ind.A]
    obs.out.B <- y[ind.B]
    obs.ben <- obs.out.A-obs.out.B
  }
  
  # Matching plot
  if (MatchingPlot== TRUE) {
    plot(pred.ben.A, pred.ben.B)
    lines(c(-1,1),c(-1,1),lwd=2)
  }
  
  # Benefit c-statistic
  cindex <- rcorr.cens(pred.ben.avg, obs.ben)
  c.benefit <- cindex["C Index"][[1]]
  c.benefit.se <- cindex["S.D."][[1]]/2 # Half the sd of Dxy
  c.benefit
  # [1] 0.6168794 # The c-for-benefit
  c.benefit - 1.96*c.benefit.se
  # [1] 0.5742719 # The 95% lower bound of the c-for-benefit
  c.benefit + 1.96*c.benefit.se
  # [1] 0.6594868 # The 95% upper bound of the c-for-benefit
  cat("c-for-benefit: ", round(c.benefit, 2), "[", round(c.benefit - 1.96*c.benefit.se,2),
      "; ",  round(c.benefit + 1.96*c.benefit.se,2), "]", sep="")
  return(list(c.benefit, c.benefit.se))
}    

find_bin <- function(vals, up, lo) {
  res <- c()
  for (i in 1:length(vals)) {
    
    for (j in 1:length(up)) {
      if (vals[i] >= lo[j] & vals[i] <= up[j]) res <- c(res, j)
    }
  }
  return(res)
}

summarise.by.group <- function(response, group, func) {
  response.split <- split(response, group)
  sapply(response.split, func)
}

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

calibration <- function(pred_rr, y, w, t, cens_time, n_bins=5) {
  # Form a calibration plot of predicted risk reduction, return the slope, intercept, predicted and observed ARRs
  
  nonNAs <- which(is.na(pred_rr)==FALSE)
  pred_rr <- pred_rr[nonNAs]
  y <- y[nonNAs]
  w <- w[nonNAs]
  t <- t[nonNAs]
  
  bins <- quantile(pred_rr, probs=linspace(0, 1, n_bins + 1))
  bins2 <- bins(as.numeric(pred_rr), n_bins, minpts=20)
  bins3 <- bins.getvals(bins2, minpt = -Inf, maxpt = Inf)

  quantiles <- find_bin(pred_rr, up=attributes(bins3)$binhi, lo=attributes(bins3)$binlo)
  #quantiles = np.digitize(pred_rr, bins) - 1
  #quantiles[quantiles == n_bins] = n_bins - 1
  tempdf <- data.frame(pred_rr=pred_rr, quantiles=quantiles, y=y, w=w, t=t)
  mean_pred <- summarise.by.group(tempdf$pred_rr, tempdf$quantiles, mean)
  tempdf$pred_rrs <- 0
  for (i in 1:nrow(tempdf)) {
    tempdf$pred_rrs[i] <- mean_pred[which(as.numeric(names(mean_pred)) == tempdf$quantiles[i])]
  }#  pred_rr = [np.mean(pred_rr[quantiles == i]) for i in range(n_bins)]

  tempdf <- tempdf[order(tempdf$quantiles),]
  
  obs_rr = c()
  for (i in 1:n_bins) {
    with_rx = which(tempdf$quantiles == i & tempdf$w == 1)
    no_rx = which(tempdf$quantiles == i & tempdf$w == 0)
    
    kmf_with_rx = with(tempdf[with_rx,], Surv(t, y))
    kmf_no_rx = with(tempdf[no_rx,], Surv(t, y))
    #kmf_no_rx = KaplanMeierFitter().fit(t[no_rx], y[no_rx])
    
    kmfit_with_rx <- survfit(Surv(t, y)~1, data=tempdf[with_rx,])#    surv = kmf_with_rx.survival_function_.reset_index()
    idx = R_searchsorted(kmfit_with_rx$time, v=cens_time) # Find indices where elements should be inserted to maintain order
    with_rx <- 1 - summary(kmfit_with_rx)$surv[min(idx, length(summary(kmfit_with_rx)$surv) - 1)]
    
    kmfit_no_rx <- survfit(Surv(t, y)~1, data=tempdf[no_rx,])#    surv = kmf_with_rx.survival_function_.reset_index()
    idx = R_searchsorted(kmfit_no_rx$time, v=cens_time) # Find indices where elements should be inserted to maintain order
    no_rx <- 1 - summary(kmfit_no_rx)$surv[min(idx, length(summary(kmfit_no_rx)$surv) - 1)]

    obs_rr <- c(obs_rr, (no_rx - with_rx))
  }
  
  #lm_df <- data.frame(pred_rr=tempdf$pred_rr, obs_rr=obs_rr)
  lm_df <- data.frame(pred_rr=mean_pred, obs_rr=obs_rr)
  cal_lm <- lm(obs_rr ~ pred_rr, data=lm_df)
  cal_slope <- as.numeric(cal_lm$coefficients[2])
  cal_intercept <- as.numeric(cal_lm$coefficients[1])
  #pred_rr, obs_rr = np.array(pred_rr), np.array(obs_rr)
  #slope, intercept = np.polyfit(pred_rr, obs_rr, 1) #Least squares polynomial fit, 1 degree of fitting
  return(list(cal_slope, cal_intercept, pred_rr, obs_rr))
}



base_dir = paste(hauptpfad, "/results/", model, "/", data.source, sep="")
setwd(base_dir)

pred_rr <- readRDS(paste0(base_dir, "/pred_rr.rds")); pred_rr_R <- pred_rr
X <- readRDS(paste0(base_dir, "/X.rds"))
w_R <- readRDS(paste0(base_dir, "/w.rds")); w <- w_R
y <- readRDS(paste0(base_dir, "/y.rds"))
t <- readRDS(paste0(base_dir, "/t.rds"))
y_cut <- readRDS(paste0(base_dir, "/y_cut.rds"))
cens <- readRDS(paste0(base_dir, "/cens.rds"))
cens_ind <- which(all_data[["cens"]] == 0)
#dataset_cut = data.frame(X= X[cens == 0,], w=w[cens == 0], y=y_cut[cens == 0])

# read python results
  if (python==TRUE) {
    require(RcppCNPy)
    pred_rr_python <- npyLoad(paste0(base_dir, "/pred_rr.npy")) # pred_rr <- pred_rr_python
    y_python <- npyLoad(paste0(base_dir, "/y.npy"), type="integer")
    w_python <- npyLoad(paste0(base_dir, "/w.npy"))
    y_cut_python <- npyLoad(paste0(base_dir, "/y_cut.npy"))
  }

# similar results, if 
  # pred_rr <- pred_rr * (-1)
  # ???

# metrics for the dataset, evaluated as dichotomous outcome
  arr_ben <- arr_noben <- c_stat <- c()
  set.seed <- 123
  print("Calculate c statistic for benefit in bootstrap samples")
  pb <- txtProgressBar(min=0, max=bootstrap_samples, style=3)
  
  
  evaluate_bootstrap_disc <- function(metric="c_for_benefit") {
    idxs = 1:nrow(dataset[["X"]])
    treated_cvd = which(dataset[["w"]] == 1 & dataset[["y"]] == 1)
    treated_nocvd = which(dataset[["w"]] == 1 & dataset[["y"]] == 0)
    untreated_cvd = which(dataset[["w"]] == 0 & dataset[["y"]] == 1)
    untreated_nocvd = which(dataset[["w"]] == 0 & dataset[["y"]] == 0)
    idxs = c(sample(treated_cvd, size=length(treated_cvd), replace=TRUE),
             sample(treated_nocvd, size=length(treated_nocvd), replace=TRUE),
             sample(untreated_cvd, size=length(untreated_cvd), replace=TRUE),
             sample(untreated_nocvd, size=length(untreated_nocvd),replace=TRUE))
    
    cens_pred_rr <- pred_rr[cens_ind]
    cens_y_cut <- y_cut[cens_ind]
    cens_w <- w[cens_ind]
    
    bucketed_arr <- bucket_arr(cens_pred_rr[idxs],cens_y_cut[idxs], cens_w[idxs])

    if (metric=="c_for_benefit") {
      ret <- c_statistic(cens_pred_rr[idxs],cens_y_cut[idxs],cens_w[idxs])
    } else if (metric=="arr_ben") {
      ret <- as.numeric(bucketed_arr[[1]])
    } else if (metric=="arr_noben") { 
      ret <- as.numeric(bucketed_arr[[2]])
    } else if (metric=="arr_noben") {
      
    }
    return(ret)
  }
  
  for (b in 1:bootstrap_samples) {
    
    setTxtProgressBar(pb, value=b)
    temp <- bucketed_arr_ben()
    arr_ben <- c(arr_ben, bucketed_arr_ben())
    
  }
  
  rmst <- slope <- intercept <- c()
  
  c_stat <- unlist(replicate(250, evaluate_bootstrap_disc("c_for_benefit"), simplify=FALSE))
  arr_ben <- unlist(replicate(250, evaluate_bootstrap_disc("arr_ben"), simplify=FALSE))
  arr_noben <- unlist(replicate(250, evaluate_bootstrap_disc("arr_noben"), simplify=FALSE))
  
  for (b in 1:bootstrap_samples) {
    
    setTxtProgressBar(pb, value=b)
    
    #idxs<- t(mapply(b, function(x) bootstrap_dataset(dataset_cut)))
    
    #idxs = bootstrap_dataset(dataset_cut); #print(head(idxs))
    idxs = 1:nrow(dataset[["X"]])
    treated_cvd = which(dataset[["w"]] == 1 & dataset[["y"]] == 1)
    treated_nocvd = which(dataset[["w"]] == 1 & dataset[["y"]] == 0)
    untreated_cvd = which(dataset[["w"]] == 0 & dataset[["y"]] == 1)
    untreated_nocvd = which(dataset[["w"]] == 0 & dataset[["y"]] == 0)
    idxs = c(sample(treated_cvd, size=length(treated_cvd), replace=TRUE),
             sample(treated_nocvd, size=length(treated_nocvd), replace=TRUE),
             sample(untreated_cvd, size=length(untreated_cvd), replace=TRUE),
             sample(untreated_nocvd, size=length(untreated_nocvd),replace=TRUE))

    cens_pred_rr <- pred_rr[cens_ind]
    cens_y_cut <- y_cut[cens_ind]
    cens_w <- w[cens_ind]
    
    bucketed_arr <- bucket_arr(cens_pred_rr[idxs],cens_y_cut[idxs], cens_w[idxs])
    #print(bucketed_arr[[1]])
    arr_ben <- c(arr_ben, bucketed_arr[[1]] )
    arr_noben <- c(arr_noben, bucketed_arr[[2]] )
    
    
    
    c_stat <- c(c_stat, c_statistic(cens_pred_rr[idxs],cens_y_cut[idxs],cens_w[idxs]))

  }
  print(get_range(c_stat)); print(get_range(arr_ben)); print(get_range(arr_noben))
  # python results
  # arr_ben: [0.03, 0.03, 0.04]#  arr_noben: [-0.04, -0.03, -0.01]
  
  
# metrics for the dataset, evaluated on the entire sample
  rmst <- slope <- intercept <- c()
  set.seed <- 123
  print("Calculate calibration in bootstrap samples")
  pb <- txtProgressBar(min=0, max=bootstrap_samples, style=3)
  for (b in 1:bootstrap_samples) {
    setTxtProgressBar(pb, value=b)
    
    idxs = bootstrap_dataset(dataset_all)
    
    rmst <- c(rmst, decision_value_rmst(pred_rr[idxs], y[idxs],
                                        w[idxs], t[idxs],
                                        args.cens_time))

    cali <- calibration(pred_rr[idxs], y[idxs],w[idxs], t[idxs],args.cens_time, n_bins=5)
    slope <- c(slope, cali[[1]]); intercept <- c(intercept, cali[[2]])
  }
  print(get_range(rmst)); print(get_range(slope))#; print(get_range(intercept));

# metrics for the dataset, non-bootstrapped
  pvals = wald_test(pred_rr[cens == 0], y_cut[cens == 0], w[cens == 0])

  stats <- list(c_stat=c_stat, rmst=rmst, slope=slope, intercept=intercept,
                arr_ben=bucketed_arr[[1]],arr_noben=bucketed_arr[[2]],
                p_ben=pvals[[1]], p_noben=pvals[[2]])
  
  for (i in 1:length(stats)) {
    datname <- names(stats)[i]
    ob <- stats[[i]]
    pfad <- paste0(base_dir, "/", datname, ".rds")
    saveRDS(ob, pfad)
    #save(ob, file=pfad)
    #eval(parse(text=paste("save(ob, file=paste0(base_dir, '/", datname, ".RData'))",sep="")))
  }
  cat("== Saved to: ", base_dir, sep="")
  

