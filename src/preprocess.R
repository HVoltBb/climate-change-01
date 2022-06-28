"
Author/maintainer: Can Zhou [eidotog@gmail.com]
Date: May 29 2021
Version: 0.3
"

require(mgcv)

# Use the hash code below to check for changes/updates in the database
# dbe6ad1720b5c377b57bc090281ab561
# data_version_hash = tools::md5sum('Replace_me_with_your_downloaded_database_file_name.csv')

lag_= 0

dat.c = dat[,c('ANIMALID', 'REL_VESS_NATION', 'REL_YEAR', 'REL_MONTH', 'REL_DAY', 'REL_SPECIES_ITIS', 'REL_SEX', 'REC_SEX','REL_LENGTH', 'REC_YEAR', 'REC_MONTH', 'REC_DAY', 'REC_LENGTH')]
dat.c = dat.c[!is.na(dat.c$REL_LENGTH) & !is.na(dat.c$REC_LENGTH) & !is.na(dat.c$REL_YEAR) & !is.na(dat.c$REC_YEAR)& !is.na(dat.c$REL_MONTH) & !is.na(dat.c$REC_MONTH),]
dat.c$ANIMALID = as.numeric(substr(dat.c$ANIMALID, 1, 12))

spp_names = unique(dat.c$REL_SPECIES_ITIS)
dat.c$spp = 0
for(i in 1:length(spp_names)){
    dat.c$spp[which(dat.c$REL_SPECIES_ITIS==spp_names[i])] = i   
}

# Date time conversion to real numbers ----
t1 = as.Date(paste(dat.c$REL_YEAR, dat.c$REL_MONTH, dat.c$REL_DAY, sep='-'), "%Y-%m-%d") - lag_
ym = zoo::as.yearmon(t1)
t1_next_month = zoo::as.Date(ym, frac=1) + 1
t2 = as.Date(paste(dat.c$REC_YEAR, dat.c$REC_MONTH, dat.c$REC_DAY, sep='-'), "%Y-%m-%d") - lag_
ym2 = zoo::as.yearmon(t2)
t2_next_month = zoo::as.Date(ym2, frac=1) + 1

d1 = as.numeric(difftime(t1_next_month, t1, units = 'days'))
d2 = as.numeric(difftime(t2_next_month, t2, units = 'days'))

dat.c$d1 = d1/30
dat.c$d2 = d2/30

# Sexing errors ----
# SEX code: 0, unknown; 1, male; 2, female
dat.c$SEX = dat.c$REL_SEX
dx = which(dat.c$REL_SEX!=dat.c$REC_SEX)
dat.c$SEX[dx] = max(dat.c$REC_SEX[dx], dat.c$REL_SEX[dx])
sexing_errors = which(dat.c$REC_SEX!=dat.c$REL_SEX & dat.c$REC_SEX!=0 & dat.c$REL_SEX!=0)
dat.c$SEX[sexing_errors] = 0

# Climate indices ----
# Daily
nao.d$Date = 1:(dim(nao.d)[1]) # 1950-01-01 => Day 1

day1 = as.Date("1950-1-1", "%Y-%m-%d")
d1.d = as.integer(difftime(t1, day1, units = "days"))
d2.d = as.integer(difftime(t2, day1, units = "days"))

dat.c$d1.d = d1.d
dat.c$d2.d = d2.d

# Monthly
year.min = min(dat.c$REL_YEAR)-1 # 1962
dat.c$REL_Date = (dat.c$REL_YEAR-year.min)*12 + dat.c$REL_MONTH
dat.c$REC_Date = (dat.c$REC_YEAR-year.min)*12 + dat.c$REC_MONTH
nao$Date = (nao$V1-year.min)*12+nao$V2
nao.c = nao[which(nao$Date==0):dim(nao)[1],]
# Dec 1961 -> 0

# Number of days in a month
# starting from Dec 1961
n_this_month = as.Date(paste(nao$V1, nao$V2, "01", sep = "-"), "%Y-%m-%d")
n_ym = zoo::as.yearmon(n_this_month)
n_next_month = zoo::as.Date(n_ym, frac=1) + 1
n_days = as.numeric(difftime(n_next_month, n_this_month, units = "days"))/30
n_days.c = n_days[which(nao$Date==0):dim(nao)[1]]

# Dating errors ----
ex = which(dat.c$d2.d<dat.c$d1.d)
# Same day tag and recapture
same_day = which(dat.c$d1.d==dat.c$d2.d)
# Missing values
no_day = which(is.na(dat.c$d1) | is.na(dat.c$d2))
# Too long growth invertals
too_long = which(dat.c$d2.d - dat.c$d1.d > 365*30)
exx = union(union(union(ex, same_day), no_day), too_long)
dat.c = dat.c[-exx,]
indx = which(dat.c$spp==1)

datx = dat.c[indx,]
table(datx$SEX)

# Marking mutiple recaptures ----
w = table(datx$ANIMALID)
dupes = names(which(w>1))
singl = names(which(w==1))

duplicity = rep(0, dim(datx)[1])
duplicity[-which(datx$ANIMALID %in% dupes)]= 1:length(singl)

for(i in 1:length(dupes)){
    tmp = which(datx$ANIMALID==dupes[i])
    duplicity[tmp] = i+length(singl)
}

# Cubic Spline ----
gam_design = mgcv::gam(V4~s(V4, bs='cs', k=11), data = nao.d, fit=FALSE)
gam_design$smooth[[1]]$xp
S = gam_design$smooth[[1]]$S[[1]]
Sdim = dim(S)[1]

nao.p = c(0, seq(min(nao.d$V4, na.rm = TRUE), max(nao.d$V4, na.rm=TRUE), length.out = 100))
nao.p = sort(nao.p)
P_design = mgcv::PredictMat(gam_design$smooth[[1]], data = data.frame(V4 = nao.p))
p_dm =     mgcv::PredictMat(gam_design$smooth[[1]], data = data.frame(V4 = nao$V3))
X_m =     mgcv::PredictMat(gam_design$smooth[[1]], data = data.frame(V4 = nao.c$V3))
