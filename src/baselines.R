run_for_framingham <- function(dataset) {
  risks <- rep(0, length(dataset[["y"]]))
  pb<-txtProgressBar(min=0, max=nrow (dataset[["X"]]), style=3)
  for (i in 1:nrow (dataset[["X"]])) {
    setTxtProgressBar(pb, value=i)
    gender <- ifelse(dataset[["X"]][i,"FEMALE"] == 1, "F", "M")
    age = dataset[["X"]][i,"AGE"]
    bmi = dataset[["X"]][i,"BMI"]
    sbp = dataset[["X"]][i,"SBP.y"]
    smk = dataset[["X"]][i,"currentsmoker"]
    dia = dataset[["X"]][i,"diabetes"]
    risks[i] <- frs(gender, 10, age=age, bmi=bmi, sbp=sbp, ht_treat=FALSE, smk=smk, dia=dia)
  }
  return(risks)
}
  
    
run_for_ascvd <- function(dataset) {
  risks <- rep(0, length(dataset[["y"]]))
  # coeffs
  male_white <- c(12.344, 11.853, -2.664, -7.99, 1.769, 1.764, 7.837, -1.795, 0.658)
  male_black <- c(2.469, 0.302, 0, -0.307, 0, 1.809, 0.549, 0, 0.645)
  female_white <- c(-29.799, 4.884, 13.540, -3.114, -13.578, 3.149, 1.957, 0, 7.574, -1.665, 0.661)
  female_black <- c(17.114, 0, 0.940, 0, -18.920, 4.475, 27.820, -6.087, 0.691, 0, 0.874)
  pb<-txtProgressBar(min=0, max=nrow (dataset[["X"]]), style=3)
  for (i in 1:nrow (dataset[["X"]])) {
    setTxtProgressBar(pb, value=i)
    female <- ifelse(dataset[["X"]][i,"FEMALE"] == 1, TRUE, FALSE)
    black <- ifelse(dataset[["X"]][i,"RACE_BLACK"] == 1, TRUE, FALSE)
    age <- log(dataset[["X"]][i,"AGE"])
    age_sq = age^2
    chol = log(dataset[["X"]][i,"CHR"])
    age_chol = age * chol
    hdl = log(dataset[["X"]][i,"HDL"])
    age_hdl = age * hdl
    sbp = log(dataset[["X"]][i,"SBP.y"])
    age_sbp = age * sbp
    smk = dataset[["X"]][i,"currentsmoker"]
    age_smk = age * smk
    dia = dataset[["X"]][i,"diabetes"]
    
    if (female == FALSE & black == FALSE) {
      vec = c(age, chol, age_chol, hdl, age_hdl, sbp, smk, age_smk, dia)
      risks[i] = 1 - 0.9144 ^ (exp((male_white%*%vec) - 61.18))
    } else if (female == FALSE & black == TRUE) {
      vec = c(age, chol, age_chol, hdl, age_hdl, sbp, smk, age_smk, dia)
      risks[i] = 1 - 0.8954 ^ (exp((male_black%*%vec) - 19.54))
    } else if (female == TRUE & black == FALSE) {
      vec = c(age, age_sq, chol, age_chol, hdl, age_hdl, sbp, age_sbp, smk,age_smk, dia)
      risks[i] = 1 - 0.9665 ^ (exp((female_white%*%vec) - 29.18))
    } else if (female == TRUE & black == TRUE) {
      vec = c(age, age_sq, chol, age_chol, hdl, age_hdl, sbp, age_sbp, smk, age_smk, dia)
      risks[i] = 1 - 0.9533 ^ (exp((female_black%*%vec) - 86.61))
    }
  }
  return(risks)
}
 

get_decision_value_rmst_naive <- function(dataset, cens_time) {
  pred_rr <- rep(0, length(dataset[["y"]]))
  pred_rr[which(dataset[["X"]][, "diabetes"] == 0)] <- 0.1
  pred_rr[which(dataset[["X"]][, "diabetes"] == 1)] <- -0.1
  return(decision_value_rmst(pred_rr, dataset[["y"]], dataset[["w"]], dataset[["t"]], cens_time))
}

    base_dir = paste(hauptpfad, "/results/baselines/", data.source, sep="")
    setwd(base_dir)
    if (calc_baseline_risk == TRUE) {
      print("Calculating baseline risks...")
      frs_res <- run_for_framingham(all_data)
      ascvd_res <- run_for_ascvd(all_data)
      #save(frs_res, file=paste0(base_dir, "/framingham.RData"))
      #save(ascvd_res, file=paste0(base_dir, "/ascvd.RData"))
      saveRDS(frs_res, paste0(base_dir, "/framingham.rds"))
      saveRDS(ascvd_res, paste0(base_dir, "/ascvd.rds"))
    }
    if (calc_naive_rmst == TRUE) {
      print("Calculating naive strategy RMST...")
      rmsts = c()
      set.seed <- 123
      pb<-txtProgressBar(min=0, max=bootstrap_samples, style=3)
      for (i in 1:bootstrap_samples) {
        setTxtProgressBar(pb, value=i)
        idxs <- bootstrap_dataset(all_data)
        alldat.X <- all_data[["X"]]; alldat.X <- alldat.X[idxs,]
        alldat.y <- all_data[["y"]]; alldat.y <- alldat.y[idxs]
        alldat.w <- all_data[["w"]]; alldat.w <- alldat.w[idxs]
        alldat.t <- all_data[["t"]]; alldat.t <- alldat.t[idxs]
        
        bootdata <- list(X=alldat.X, y=alldat.y, w=alldat.w, t=alldat.t)
        rmsts <- c(rmsts, (get_decision_value_rmst_naive(bootdata, cens_time= 365.25 * 3)))
      }
      print(get_range(rmsts))
    }
  