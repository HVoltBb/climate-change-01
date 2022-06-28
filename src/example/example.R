"
Author/maintainer: Can Zhou [eidotog@gmail.com]
Date: May 29 2021
Version: 0.1
"
require(TMB)

dat = read.csv('../../data/MRMVRELRECAP03.csv')
nao.d = read.csv("../../data/norm.daily.nao.index.b500101.current.ascii", sep="", header = FALSE)
nao = read.csv('../../data/norm.nao.monthly.b5001.current.ascii', sep='', header=FALSE)
source('../preprocess.R')

TMB::compile('example.cpp')
dyn.load(TMB::dynlib('example'))

pickone = 1

obj = TMB::MakeADFun(
    data = list(l1 = datx$REL_LENGTH,
                l2 = datx$REC_LENGTH,
                t1 = datx$REL_Date,
                t2 = datx$REC_Date,
                d1 = datx$d1,
                d2 = datx$d2,
                days = n_days.c,
                PAR1 = 0 # 0: von Bertalanffy, 1: G. West, 2: Gompertz, 3: Logistic, 4: Generalized logistic, 5+: General von Bertalanffy
    ),
    parameters  = list( c = c(250, 60)[pickone],
                        k = c(0.05, 0.06)[pickone],
                        logSig = 3,
                        lognu = 0),
    map = list(
        lognu = factor(NA),
        c = factor(1)), # Do not edit this line. Just a placeholder.
    DLL = 'example',
    silent = FALSE
)

fit=nlminb(obj$par, obj$fn, obj$gr)
# AIC
2*(fit$objective+length(fit$par)) 
# 72640.17 for von Bertalanffy growth

dyn.unload(TMB::dynlib('example'))
