"
Author/maintainer: Can Zhou [eidotog@gmail.com]
Date: May 29 2021
Version: 0.3
"
# Libs ----
require(TMB)

# Data section ----
dat <- read.csv("../data/MRMVRELRECAP03.csv")
nao.d <- read.csv("../data/norm.daily.nao.index.b500101.current.ascii", sep = "", header = FALSE)
nao <- read.csv("../data/norm.nao.monthly.b5001.current.ascii", sep = "", header = FALSE)
source("preprocess.R")

# Compile and load the program ----
TMB::compile("v5.cpp", flags = "-D_ZOOM -Wno-ignored-attributes")
dyn.load(TMB::dynlib("v5"))

# Model fitting ----
pickone <- 2 # Pick a set of inital values to start with

obj <- TMB::MakeADFun(
  data = list(
    l1 = datx$REL_LENGTH,
    l2 = datx$REC_LENGTH,
    t1 = datx$REL_Date,
    t2 = datx$REC_Date,
    d1 = datx$d1,
    d2 = datx$d2,
    d1_d = datx$d1.d,
    d2_d = datx$d2.d,
    SEX = datx$SEX,
    DUPE = duplicity - 1,
    X = X_m,
    X_d = gam_design$X[, -1],
    p_dm = p_dm,
    S = Matrix::.bdiag(list(S)),
    Sdim = dim(S)[1],
    Pdim = dim(P_design)[1],
    prediction_design_matrix = P_design,
    MNS = datx$REL_MONTH - 1,
    MNS_map = nao.d$V2 - 1,
    nao = c(nao.c$V3),
    nao_d = nao.d$V4,
    days = n_days.c,
    days_p = n_days[-1],
    lag_ = lag_,
    REFSEX = 2,
    REFLEN = 150.0,
    REFTIM = 1.0,
    REFMON = 0,
    PAR1 = 1, # 0: von Bertalanffy, 1: G. West, 2: Gompertz, 3: Logistic, 4: Generalized logistic, 5+: General von Bertalanffy
    PAR2 = 0, # 0: additive effect, 1: multiplicative effect. Deprecated since v5
    PAR3 = .0, # Only considered when PAR1 >= 5. Deprecated since v5
    PAR4 = 60, # break point between daily resolution and monthly resolution
    PAR5 = 1, # 1: linf, 0: k [location of intrinsic effects]
    PAR6 = 0 # 0: no measurement error, otherwise, with measurement error
  ),
  parameters = list(
    c = c(240, 60.6)[pickone],
    b = 0,
    k = c(0.04, 0.058)[pickone],
    logSig = 3.2,
    logSexl = 0,
    logSexk = -4.85,
    logNaoc = 0,
    logNaok = 8.78,
    logMonthl = 0,
    lognu = 0,
    logvIndl = -1,
    lsex = rep(0, 3),
    ksex = rep(0, 3),
    gc = rep(0, dim(S)[1]),
    gk = rep(0, dim(S)[1]),
    theta_l = rep(0, 12),
    e_o = rep(0, dim(datx)[1]),
    indl = rep(0, length(dupes) + length(singl))
  ),
  map = list(
    # Instructions:
    # Comment OFF the each line to turn the corresponding effect ON,
    # either in combination or individually except the month effect
    # and seasonal pattern, which at most one of them can be turned
    # ON simultaneously.
    logSexk = factor(NA), ksex = factor(rep(NA, 3)), # sex effect on k
    logSexl = factor(NA), lsex = factor(rep(NA, 3)), # sex effect on Linf

    indl = factor(rep(NA, length(dupes) + length(singl))), logvIndl = factor(NA), # individual effect
    lognu = factor(NA),
    e_o = factor(rep(NA, dim(datx)[1])), # Observational error
    logMonthl = factor(NA), theta_l = factor(rep(NA, 12)), # month effect

    b = factor(NA), # linear effect
    # logNaoc = factor(NA), gc = factor(rep(NA, dim(S)[1])), # non-linear effect on l
    logNaok = factor(NA), gk = factor(rep(NA, dim(S)[1])), # non-linear effect on k
    c = factor(1)
  ), # Do not edit this line. Just a placeholder.
  random = c("lsex", "ksex", "gc", "theta_l", "indl", "gk", "e_o"),
  DLL = "v5",
  silent = FALSE # You may want to change it to FALSE to check for errors
)

fit <- nlminb(fit$par, obj$fn, obj$gr)

# Postprocessing ----
# AIC
2 * (fit$objective + length(fit$par))
# This should give you the same result as in the example code when this code is unmodified.

# Fixed parameters and their standard errors
rep <- TMB::sdreport(obj)

# Random parameters and their standard errors
rep$par.random[names(rep$par.random) == "ksex"]
sqrt(rep$diag.cov.random[names(rep$par.random) == "ksex"])
rep$par.random[names(rep$par.random) == "gk"]

dyn.unload(TMB::dynlib("v5"))