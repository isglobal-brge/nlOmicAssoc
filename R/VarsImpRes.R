dataDir<-"/home/lnonell/nlOmicAssoc_tests"
resultsDir<-"/home/lnonell/nlOmicAssoc_tests/Results"

require(Biobase)
require(mice)
load(file.path(dataDir,"gexp4b.Rdata")) #object gexp4b
data_m_all<-exprs(gexp4b)
exps<-read.table(file=file.path(dataDir,"exposures.tsv"),sep = "\t",dec=".",stringsAsFactors = F,header = T)
exps_desc<-read.table(file=file.path(dataDir,"description.tsv"),sep = "\t",dec=".",stringsAsFactors = F,header = T)
#select pregnancy variables
vars_preg<-exps_desc[!is.na(exps_desc$TimePoint) & exps_desc$TimePoint=="Pregnancy",]

load(file.path(dataDir,"20170510_ImputedDataset.Rdata")) #object imp
source(file.path(dataDir,"R/Fit_LN_DEF.R"))

#data
data_m<-data_m_all

#aquestes són les variables amb NAs
vars_df<-exps[,-1] #trec la primera i li poso a l'id
rownames(vars_df) <- exps[,1]

#variables imputades per Ibon
vars_df<-complete(imp)[,-112] #last is id

Fit(data_m, vars_df, df = 4, select = 0.05, select_adj = NULL, cores = 1L, imp.method='fastpmm')


library(lattice)
pdf(file=file.path(resultsDir,"vars.density.imputedIbon.vs.original.pdf")) 
for (i in c(1:16,18:ncol(vars_df))){
   var<-colnames(vars_df)[i];
   print(densityplot(imp, as.formula(paste0("~",var)),main=var))
 }
dev.off()

st1<-Sys.time()
fit.all<-Fit(data_m,vars_df,cores=12)
st2<-Sys.time()
st2-st1

##################################################
#sets de proves
vars_df<-exps[,-1] #trec la primera i li poso a l'id
rownames(vars_df) <- exps[,1]
vars_df<-complete(imp)[,-112] #last is id
vars_df<-vars_df[,vars_preg$Exposure] #dades imputades
data_m<-data_m_all[1:50,]
#faig córrer la funció step by step amb els següents params
df = 4; select = 0.05; select_adj = NULL; cores = 1L; imp.method='fastpmm'
fit<-Fit(data_m,vars_df,cores=4)
head(fit$table)
head(fit$selected_vars)
fit$analysis_vars
fit$alpha_select
fit.s1<-Fit(data_m,vars_df,select_adj=1,cores=4)
head(fit$table)
head(fit.s1$selected_vars)
fit.s1$analysis_vars
fit.s1$alpha_select
#són molt diferents els resultats
fit.s05<-Fit(data_m,vars_df,select_adj=0.05,cores=4)
head(fit.s05$table)
head(fit.s05$selected_vars)
fit.s05$analysis_vars
fit.s05$alpha_select


# save(fit.all,file=file.path(resultsDir,"file.all.RData"))
# load(file=file.path(resultsDir,"file.all.RData"))
# head(fit.all$table)
# fit.all$selected_vars
# fit.all$analysis_vars
# fit.all0$alpha_select

######Methylation
load(file.path(dataDir,"methy4b.Rdata")) #object methy4b
methy4b
#Expression set 465.930 x 308 samples
meth.all<-exprs(methy4b)
hist(meth.all)
#we need to filter data before applying mfp
#study sd
meth.all.sd<-apply(meth.all,1,sd)
hist(meth_all.sd)

meth.all.f<-meth.all[meth.all.sd>0.1,]
dim(meth.all.f)
meth.all.f.sd<-apply(meth.all.f,1,sd)
