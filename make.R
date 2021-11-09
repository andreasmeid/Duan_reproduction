# set path
  laufwerk <- "J"
  hauptpfad <- paste0(laufwerk, ":/Duan/hte-prediction-rcts-master")

# load functions
  source(as.character(paste0(hauptpfad, "/src/framingham.R")))
  source(as.character(paste0(hauptpfad, "/src/dataloader.R")))

# set arguments
  args.dataset <- data.source <- "combined"         # which trial data to use
  calc_baseline_risk <- TRUE                        # if baseline risk is to be calculated
  args.calc_naive_rmst <- calc_naive_rmst <- TRUE   # if naive (baseline) rmst is to be calculated
  bootstrap_samples <- 250                          # number of bootstrap samples to draw
  args.model <- "xlearner"                          # type of predictive algorithm
  cens_time <- args.cens_time <-  365.25 * 3        # time after which patients are censored
  load_new <- FALSE                                 # if data is to be processed from the beginning
  python <- FALSE                                   # if python results are to be loaded
  args.cstat <- cstat <- TRUE                       # if c for benefit is to be calculated
  args.calc_rmst <- calc_rmst <- TRUE               # if rmst is to be calculated
  args.calibrat <- calibrat <- TRUE                 # if calibration is to be calculated
  
# load data
  if (load_new) {
    cat("== Running for: ", args.dataset, "\n", sep="")
    setwd(paste0(laufwerk, ":/Duan/hte-prediction-rcts-master/src"))
    if (args.dataset == "combined") {
      dataset <- combine_datasets(load_data("sprint"), load_data("accord"))
    } else {
      dataset <- load_data(args.dataset)
    }
    temp <- cut_dataset(dataset, cens_time=args.cens_time)
    cut_data <- dataset_cut <- bin_data <- temp[[1]]; all_data <- dataset_all <- temp[[2]]
    save.image("make_dataset.RData")
    load_new <- FALSE
  } else {
    setwd(paste0(laufwerk, ":/Duan/hte-prediction-rcts-master/src"))
    load("make_dataset.RData"); 
    laufwerk <- "J"
    hauptpfad <- paste0(laufwerk, ":/Duan/hte-prediction-rcts-master")
  }

# calculate baselines 
  source(as.character(paste0(hauptpfad, "/src/baselines.R")))

# predict with xlearner
  #icpw <- read.table(paste0(hauptpfad, "/icpw.csv"), header=FALSE, sep="\n")$V1
  #cut_data[["ipcw"]] <- icpw[cut_data[["X"]][,1]]
  source(paste0(hauptpfad, "/lib/xlearner-rf.R"))
  args.model <- model <- "xlearner"
  source(paste0(hauptpfad, "/src/predict.R"))
# predict with logistic regression     
  args.model <- model <- "logreg"
  source(paste0(hauptpfad, "/src/predict.R"))
  
# evaluate xlearner 
  #python <- TRUE
  #bootstrap_samples <- 250 
  args.model <- model <- "xlearner"
  source(paste0(hauptpfad, "/src/evaluate.R"))
  
