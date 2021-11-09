# Calculates difference between performance on a bootstrapped dataset (upon which the model is trained) 
# and the original dataset. Optimism is defined as the mean difference over many bootstrap datasets.



run_for_optimism <- function (original_dataset, bootstrap_dataset, args.model="logreg") {
  boot_cutted <- cut_dataset(bootstrap_dataset, args.cens_time)
  cut_data <- boot_cutted[[1]]; all_data <- boot_cutted[[2]] 
  orig_cutted <- cut_dataset(original_dataset, args.cens_time)
  cut_data_orig <- orig_cutted[[1]]; all_data_orig <- orig_cutted[[2]]
  
  if(args.model == "cox") {
    model = CoxAIC()
    model.train(all_data["X"], all_data["w"], all_data["y"], all_data["t"])
    pred_rr = model.predict(args.cens_time, all_data["X"])
    pred_rr_orig = model.predict(args.cens_time, all_data_orig["X"])
  } else if (args.model == "logreg") {
    X <- get_interaction_terms(X=cut_data[["X"]], w=cut_data[["w"]])
    X <- cbind(X, cut_data[["w"]]); names(X)[ncol(X)] <- "w"
    y <- cut_data[["y"]]; ipcw <- cut_data[["ipcw"]]
    train.logreg <- glm(data = X, y ~ ., weights=ipcw)
    
    X_all = get_interaction_terms(X=all_data[["X"]], w=rep(1, nrow(all_data[["X"]])))
    X_all = cbind(X_all, rep(1, nrow(all_data[["X"]]))); names(X_all)[ncol(X_all)] <- "w"
    pred_rr = predict(train.logreg, newdata = X_all, type='response')

    X_orig = get_interaction_terms(X=all_data_orig[["X"]], w=rep(1, nrow(all_data_orig[["X"]])))
    X_orig = cbind(X_orig, rep(1, nrow(all_data_orig[["X"]]))); names(X_orig)[ncol(X_orig)] <- "w"
    pred_rr_orig = predict(train.logreg, newdata = X_orig, type='response')
  } else if (args.model == "linearxlearner") {
    model = LinearXLearner()
    model.train(cut_data["X"], cut_data["w"], cut_data["y"],
                cut_data["ipcw"])
    pred_rr = model.predict(all_data["X"], all_data["w"], False)
    pred_rr_orig = model.predict(all_data_orig["X"],
                                 all_data_orig["w"], False)
  }
  cens_ind <- which(cut_data[["cens"]] == 0)
  cens_pred_rr <- pred_rr[cens_ind]
  cens_y_cut <- cut_data[["y_cut"]][cens_ind]
  cens_w <- cut_data[["w"]][cens_ind]
  
  cens_ind_orig <- which(cut_data_orig[["cens"]] == 0)
  cens_pred_rr_orig <- pred_rr_orig[cens_ind_orig]
  cens_y_cut_orig <- cut_data_orig[["y_cut"]][cens_ind_orig]
  cens_w_orig<- cut_data_orig[["w"]][cens_ind_orig]
  
  c_stat_bootstrap = c_statistic(cens_pred_rr, cens_y_cut, cens_w)
  c_stat_original = c_statistic(pred_rr_orig,  cens_y_cut_orig, cens_w_orig)
  
  #c_stat <- c(c_stat, c_statistic(cens_pred_rr[idxs],cens_y_cut[idxs],cens_w[idxs]))
  
  rmst_bootstrap = decision_value_rmst(pred_rr, all_data["y"], all_data["w"], all_data["t"], args.cens_time)
  rmst_original = decision_value_rmst(pred_rr_orig, all_data_orig["y"],
                                      all_data_orig["w"], all_data_orig["t"],
                                      args.cens_time)
  
  return(list(c_stat_diff=c_stat_bootstrap - c_stat_original,
              decision_value_rmst_diff=rmst_bootstrap - rmst_original))
}


    c_stat_diff <- decision_value_rmst_diff <- c()

    evaluate_optimisn <- function(typus="c_bene") {
      
      idxs = 1:nrow(dataset_all[["X"]])
      treated_cvd = which(dataset_all[["w"]] == 1 & dataset_all[["y"]] == 1)
      treated_nocvd = which(dataset_all[["w"]] == 1 & dataset_all[["y"]] == 0)
      untreated_cvd = which(dataset_all[["w"]] == 0 & dataset_all[["y"]] == 1)
      untreated_nocvd = which(dataset_all[["w"]] == 0 & dataset_all[["y"]] == 0)
      idxs = c(sample(treated_cvd, size=length(treated_cvd), replace=TRUE),
               sample(treated_nocvd, size=length(treated_nocvd), replace=TRUE),
               sample(untreated_cvd, size=length(untreated_cvd), replace=TRUE),
               sample(untreated_nocvd, size=length(untreated_nocvd),replace=TRUE))
      
      bootstrap <- list(X=dataset[["X"]][idxs,],
                        y=dataset[["y"]][idxs],
                        w=dataset[["w"]][idxs],
                        t=dataset[["t"]][idxs])
      if (typus=="c_bene") {
        ret <- run_for_optimism(original_dataset=dataset, bootstrap_dataset=bootstrap)[[1]]
      } else {
        ret <- run_for_optimism(original_dataset=dataset, bootstrap_dataset=bootstrap)[[2]]
      }
      return(ret)
    }
    
    #set.seed <- 123
    #for (b in 1:bootstrap_samples) {
    #  idxs = bootstrap_dataset(dataset_all)
    #  bootstrap <- list(X=dataset[["X"]][idxs,],
    #                    y=dataset[["y"]][idxs],
    #                    w=dataset[["w"]][idxs],
    #                    t=dataset[["t"]][idxs])
      
    #  optimismus <- run_for_optimism(original_dataset=dataset, bootstrap_dataset=bootstrap)
    #  c_stat_diff <- c(c_stat_diff,optimismus[[1]])
    #  decision_value_rmst_diff <- c(decision_value_rmst_diff,optimismus[[2]])
    #}
    
    pboptions(type = "txt", style = 3, char = "=")
    c_stat_diff <- unlist(pbreplicate(bootstrap_samples, evaluate_optimisn("c_bene")))
    
    pboptions(type = "txt", style = 3, char = "=")
    decision_value_rmst_diff <- unlist(pbreplicate(bootstrap_samples, evaluate_optimisn("rmst_diff")))
    
    
    cat("mean optimism c-statistic ", mean(c_stat_diff), "\n",
        "mean optimism rmst ", mean(decision_value_rmst_diff), "\n", sep="")
