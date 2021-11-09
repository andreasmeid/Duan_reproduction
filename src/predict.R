# fuer LogReg
  get_interaction_terms <- function(X, w) {
    df2 <- as.data.frame(X[,-1])
    X <- X[,-1]
    names(df2) <- paste(names(X), "w",  sep="_")
    for (i in 1:nrow(df2)) {
      for (j in 1:ncol(df2)) {
        df2[i,j] <- w[i] * X[i,j]
      }
    }
    return(cbind(X,df2))
  }

run_with_model <- function(dataset, cut_data, all_data,
                           args.cens_time=cens_time, args.model=model) {
  #Running a particular choice of model, saves to /results folder.
  if (is.null(cut_data)) {
    temp <- cut_dataset(dataset, cens_time=args.cens_time)
    cut_data <- temp[[1]]; all_data <- temp[[2]]
  }
  if (args.model == "xlearner") {
    train.RFXLearner <- train_x_learner(X=cut_data[["X"]],#[,-which(names(cut_data[["X"]]) == "X")], 
                                        W=cut_data[["w"]], Y=cut_data[["y"]], ipcw=cut_data[["ipcw"]],
                                        pfad=paste0(hauptpfad, "/ckpts"))
    attach(train.RFXLearner)
    pred_rr <- predict_x_learner(all_data[["X"]], all_data[["w"]], FALSE, TRUE)
    pred_rr <- pred_rr * (-1)
    detach(train.RFXLearner)
  } else if (args.model == "cox") {
    modell = CoxAIC()
    modell.train(all_data[["X"]], all_data[["w"]], all_data[["y"]], all_data[["t"]])
    pred_rr = modell.predict(args.cens_time, all_data[["X"]])
  } else if (args.model == "survrf") {
    modell = SurvRF()
    modell.train(all_data[["X"]], all_data[["w"]], all_data[["y"]], all_data[["t"]])
    pred_rr = modell.predict(args.cens_time)
  #} else if (args.model == "causalforest") {
  #  modell = CausalForest()
  #  modell.train(cut_data["X"], cut_data["w"], cut_data["y"])
  #  pred_rr = np.r_[modell.predict(),
  #                  modell.predict(all_data["X"][all_data["cens"] == 1])]
  } else if (args.model == "logreg") {
    X <- get_interaction_terms(X=cut_data[["X"]][,-which(names(cut_data[["X"]]) == "X")], w=cut_data[["w"]])
    X <- cbind(X, cut_data[["w"]]); names(X)[ncol(X)] <- "w"
    y <- cut_data[["y"]]; ipcw <- cut_data[["ipcw"]]
    train.logreg <- glm(data = X, y ~ ., weights=ipcw)

    X_1 = get_interaction_terms(X=all_data[["X"]][,-which(names(cut_data[["X"]]) == "X")], w=rep(1, nrow(all_data[["X"]])))
    X_1 = cbind(X_1, rep(1, nrow(all_data[["X"]]))); names(X_1)[ncol(X_1)] <- "w"
    py1 = predict(train.logreg, newdata = X_1, type='response')
    
    X_0 = get_interaction_terms(X=all_data[["X"]][,-which(names(cut_data[["X"]]) == "X")], w=rep(0, nrow(all_data[["X"]])))
    X_0 = cbind(X_0, rep(0, nrow(all_data[["X"]]))); names(X_0)[ncol(X_0)] <- "w"
    py0 = predict(train.logreg, newdata = X_0, type='response')
    pred_rr <- py0 - py1
    
  } else if (args.model == "linearxlearner") {
    modell = LinearXLearner()
    modell.train(cut_data[["X"]], cut_data[["w"]], cut_data[["y"]],cut_data[["ipcw"]])
    pred_rr = modell.predict(all_data[["X"]], all_data[["w"]], FALSE)
  } else {
    stop("Not a supported model.")
  }
  return(list(pred_rr=pred_rr, X=all_data[["X"]], w=all_data[["w"]], y=all_data[["y"]], t=all_data[["t"]],
              y_cut=all_data[["y_cut"]], cens=all_data[["cens"]]))
}

#cat("== Running for: ", args.dataset, sep="")
#setwd("D:/Duan/hte-prediction-rcts-master/src")
#if (args.dataset == "combined") {
#  dataset <- combine_datasets(load_data("sprint"), load_data("accord"))
#} else {
#  dataset <- load_data(args.dataset)
#}


if (is.null(dataset_cut)) {
  stats <- run_with_model(dataset, args.cens_time=cens_time, args.model=model)
  
} else {
  stats <- run_with_model(dataset, cut_data=dataset_cut, all_data=dataset_all, args.cens_time=cens_time, args.model=model)
}

base_dir = paste(hauptpfad, "/results/", model, "/", data.source, sep="")
setwd(base_dir)

for (i in 1:length(stats)) {
  datname <- names(stats)[i]
  ob <- stats[[i]]
  pfad <- paste0(base_dir, "/", datname, ".rds")
  saveRDS(ob, pfad)
  #save(ob, file=pfad)
  #eval(parse(text=paste("save(ob, file=paste0(base_dir, '/", datname, ".RData'))",sep="")))
}
cat("== Saved to: ", base_dir, sep="")
