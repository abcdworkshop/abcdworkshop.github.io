##########################################
##########################################
## Code for: Mixed Models Workshop      ##
## Wes Thompson (wes.stat@gmail.com)    ##
## July 6, 2021                         ##
##########################################
##########################################

rm(list=ls())
gc()

####################
####################
## Load libraries ##
####################
####################

library(gamm4)
library(nlme)
library(Hmisc)
library(psych)
library(mice)
library(cvTools)
library(modEvA)
library(geepack)
library(MuMIn)
library(ggplot2)
library(data.table)

####################################
####################################
## Load and manipulate nda data   ##
####################################
####################################

## Read Rds file from DEAP (ABCD NDA version 3.0 release)
#nda = readRDS("nda3.0_non_image.Rds")	
nda = readRDS("DEAP-data-download.rds")
names(nda) 

#########################################
#########################################
## Select & process data for analyses  ##
#########################################
#########################################

## Select sociodemographics
ind_demog = c(which(names(nda)=="age"),which(names(nda)=="female"), which(names(nda)=="high.educ"),which(names(nda)=="married"),which(names(nda)=="household.income"))
names(nda)[ind_demog]
summary(nda[,ind_demog])

## Select nesting variables
ind_nest = c(which(names(nda)=="rel_family_id"));summary(nda[,ind_nest])

## Select nc measures
ind_nc = c(which(names(nda)=="nihtbx_fluidcomp_uncorrected"),which(names(nda)=="nihtbx_cryst_uncorrected")); names(nda)[ind_nc]

## Select cbcl variables
ind_cbcl = c(which(names(nda)=="cbcl_scr_syn_external_t"),which(names(nda)=="cbcl_scr_syn_internal_nm")); names(nda)[ind_cbcl]

## Subset data
data = nda[,c(which(names(nda)=="src_subject_id"),which(names(nda)=="event_name"),
	ind_nest,ind_demog,which(names(nda)=="rel_relationship"),which(names(nda)=="rel_group_id"),
              ind_nc,ind_cbcl),]
data = data[data$event_name == "baseline_year_1_arm_1" | data$event_name == "2_year_follow_up_y_arm_1",]
names(data);dim(data)
table(data$event_name)

######################
######################
## A baseline model ##
######################
######################

data1 = data[data$event_name == "baseline_year_1_arm_1",]
data1 = data1[complete.cases(data1),]
dim(data1)
data1$y = scale(data1$nihtbx_fluidcomp_uncorrected)
hist(data1$y)
summary(data1$y)

lme0 = lme(y~1, random = ~1 | rel_family_id, data = data1, na.action="na.omit", method = "ML")
summary(lme0)

lme1 = lme(y~ age + female + high.educ + married, random = ~1 | rel_family_id, data = data1, na.action="na.omit", method = "ML") 
summary(lme1)

anova(lme0,lme1)
r.squaredLR(lme1,lme0)

# implied covariances & correlations
data1$fam_size = 1
data1$num = 1
for(j in unique(data1$rel_family_id)){
	data1$fam_size[data1$rel_family_id == j] = length(data1$fam_size[data1$rel_family_id == j])
	data1$num[data1$rel_family_id == j] = 1:length(data1$fam_size[data1$rel_family_id == j])
}
data2 = data1[data1$fam_size==2,]
data2 = data2[order(data2$rel_family_id,data2$num),]
dim(data2)

lme0.2 = lme(y~1, random = ~1 | rel_family_id, data = data2, na.action="na.omit")
summary(lme0.2)

VarCorr(lme0.2)

var(data2$y)
0.3724939 + 0.5754415

cov(data2$y[data2$num==1],data2$y[data2$num==2])
0.3724939

cor(data2$y[data2$num==1],data2$y[data2$num==2])
0.3724939 / var(data2$y)

######################
######################
## Nonlinear Models ##
######################
######################

gamm1 = gamm(y~ age, random = list(rel_family_id=~1), data = data1[1:1000,], na.action="na.omit") 
summary(gamm1$gam)

gamm2 = gamm(y~ poly(age,2), random = list(rel_family_id=~1), data = data1[1:1000,], na.action="na.omit") 
summary(gamm2$gam)

gamms = gamm(y~ s(age), random = list(rel_family_id=~1), data = data1[1:1000,], na.action="na.omit") 
summary(gamms$gam)
plot(gamms$gam)

anova(gamm1$lme,gamm2$lme,gamms$lme)

#################################
#################################
## Imputing Missing Covariates ##
#################################
#################################

data1 = data[data$event_name == "baseline_year_1_arm_1",]
data1$id = 1:dim(data1)[1]

# Number of multiple imputed datasets & maximum number of iterations 
n.imp = 5
n.iter = 5

var.ls <- c("age","female","high.educ","married","household.income")
dat0 = data.table(data1)
dat0 <- dat0[, var.ls, with = FALSE ]
Hmisc::describe(dat0)

# Initialize model
ini = mice( dat0, m = 1, maxit = 0 )
meth = ini$meth
meth["age"] = "norm.predict"
meth["female"] = "polyreg"
meth["high.educ"] = "polyreg"
meth["married"] = "polyreg"
meth["household.income"] = "polyreg"
pred = ini$pred

set.seed(314)
post = mice( dat0, meth = meth, pred = pred, seed = 111,
              m = 1, maxit = 0)$post
dat.imp = mice( dat0, meth = meth, pred = pred, post = post,
                 seed = 1111, m = n.imp, maxit = n.iter)
rm(dat0)
# Convert imputed data to long format
dat.mi <- complete(dat.imp, action = "long", include = TRUE)
rm(dat.imp)

names(dat.mi)[1:2] <- c("imp", "id")

# Adding variables
dat = data1[,c(1,3,9:15)] # select all unimputed variables
dat.mi <- merge( dat.mi, dat, by = "id", all.x = TRUE )
dat.mi$id = NULL
names(dat.mi)

##########################################################
## Association of factor scores with CBCL Externalizing ##
##########################################################

d_rsq = array(NA, dim = c(n.imp, 2))
gamm.results = list()
for(j in 1:n.imp){
	dat.j = dat.mi[dat.mi$imp==j,]
	ind_y = c(which(names(dat.j)=="cbcl_scr_syn_external_t")); names(dat.j)[ind_y]
	ind_demog = 2:6
	ind_nest = c(8,9)
	dat.j[,ind_y] = dat.j[,ind_y]+.1
	## Run GAMMs
	p.values = array(NA, dim=c(1,length(ind_y)))
	t.values = array(NA, dim=c(1,length(ind_y)))
	results = list()
	for(i in 1:length(ind_y)){
		show(c(j,i))
		results[[i]] = list()
		form0 = paste(names(dat.j)[ind_y[i]],"~ ")
		for(k in 1:length(ind_demog)){
			form0 = paste(form0,"+",names(dat.mi)[ind_demog[k]])
		}
		form0 = formula(form0)	
		form = paste(names(dat.j)[ind_y[i]],"~ nihtbx_fluidcomp_uncorrected")
		for(k in 1:length(ind_demog)){
			form = paste(form,"+",names(dat.mi)[ind_demog[k]])
		}
		form = formula(form)	
		ran = ~(1|rel_family_id)
		ran = formula(ran)
		gamm0 = gamm4(formula = form0, random = ran, data = dat.j[1:500,], family = Gamma(link = "log"))
		gamm1 = gamm4(formula = form, random = ran, data = dat.j[1:500,], family = Gamma(link = "log"))
		results[[i]][[1]] = gamm0
		results[[i]][[2]] = gamm1
		p.values[1,i] = summary(results[[i]][[2]]$gam)$p.pv[2]
		t.values[1,i] = summary(results[[i]][[2]]$gam)$p.t[2]
		print(summary(results[[i]][[1]]$gam))
		print(summary(results[[i]][[2]]$gam))
		print(anova(results[[i]][[2]]$gam))
		print("######################################")
	}
	d_rsq[j,] = rep(0,length(ind_y))
	for(p in 1:length(ind_y)){
		d_rsq[j,p] = summary(results[[p]][[2]]$gam)$r.sq - summary(results[[p]][[1]]$gam)$r.sq
	}
	gamm.results[[j]] = results
}	
d_rsq = apply(d_rsq,2,mean)

## Create table
cbcl.gamm.tab = as.data.frame(array(0,dim=c(12,length(ind_y)*3)))
rownames(cbcl.gamm.tab) = c("Variable",names(summary(results[[1]][[2]]$gam)$p.coeff))
names(cbcl.gamm.tab) = c(names(data1)[ind_y][1],names(data1)[ind_y][1],names(data1)[ind_y][1])
cbcl.gamm.tab[1,] = rep(c("coef   ","se    ","p_value"),length(ind_y))
ind = 0
for(k in 1:length(ind_y)){
	ind = ind + 1
	coeff = summary(gamm.results[[1]][[k]][[2]]$gam)$p.table[,1]
	for(b in 2:n.imp){
		coeff = 	coeff + summary(gamm.results[[b]][[k]][[2]]$gam)$p.table[,1]
	}
	coeff = coeff/n.imp	
	cbcl.gamm.tab[2:(12),ind] = coeff
	ind = ind + 1
	se = summary(gamm.results[[1]][[k]][[2]]$gam)$p.table[,2]
	for(b in 2:n.imp){
		se = se + summary(gamm.results[[b]][[k]][[2]]$gam)$p.table[,2]
	}
	se = se/n.imp
	imp.err = 0*se
	for(b in 1:n.imp){
		imp.err = imp.err + (summary(gamm.results[[b]][[k]][[2]]$gam)$p.table[,1]-coeff)^2 # Rubin's formula
	}
	imp.err = imp.err/n.imp	
	se = se + imp.err
	cbcl.gamm.tab[2:(12),ind] = round(se,2)
	ind = ind + 1
	cbcl.gamm.tab[2:(12),ind] = 2*(1-pnorm(abs(coeff),0,se))
}
cbcl.gamm.tab

#########################
#########################
## Longitudinal Models ##
#########################
#########################

data.l = data[1:1000,]

data.l$y = scale(data.l$nihtbx_cryst_uncorrected)
hist(data.l$y)
summary(data.l$y)

lme1 = lme(y~ age, random = ~1 | rel_family_id/src_subject_id, data = data.l, na.action="na.omit") 
summary(lme1)

p <- ggplot(data = data.l, aes(x = age, y = y, group = src_subject_id))
p + geom_line()

############################
############################
## Longitudinal Stability ##
############################
############################

data.l = data[order(data$src_subject_id,data$age),]
data.l = data.l

# select neurocog variables
ind_nc = grep("nihtbx",names(data.l))
summary(data.l[,ind_nc])

data.l$visit = 1; data.l$visit[data.l$event_name == '2_year_follow_up_y_arm_1'] = 2
id = unique(data.l$src_subject_id)
n = length(id)
data.l$age0 = NA
data.l$age_d = NA
for(i in 1:n){
	show(i)
	data.l_i = data.l[data.l$src_subject_id == id[i],]
	data.l_i$age0 = data.l_i$age[data.l_i$event_name == 'baseline_year_1_arm_1']
	data.l_i$age_d = data.l_i$age - data.l_i$age0
	data.l[data.l$src_subject_id == id[i],] = data.l_i
}

data.l$y = scale(data.l$nihtbx_cryst_uncorrected)

y.lme = lme(y~age0*age_d + visit, random=~1|rel_family_id/src_subject_id, data = data.l, na.action = "na.omit")
summary(y.lme)

VarCorr(y.lme)

long_stab = (as.numeric(VarCorr(y.lme)[2])+as.numeric(VarCorr(y.lme)[4]))/(as.numeric(VarCorr(y.lme)[2])+
		as.numeric(VarCorr(y.lme)[4])+as.numeric(VarCorr(y.lme)[5]))
show(long_stab)
