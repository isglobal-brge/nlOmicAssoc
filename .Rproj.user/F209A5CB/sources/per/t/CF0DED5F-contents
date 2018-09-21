# =============================================================================
# USAGE EXAMPLE
# =============================================================================

# Preambles
# =============================================================================

rm(list = ls())
gc(reset = TRUE)
require(Biobase)

# Load expressions an espositions:
# -----------------------------------

load("objects/expo_INMA.RData")
load("objects/expset_INMA.RData")


expo_INMA[expo_INMA<0] <- 0
 
expoINMA.log <- data.frame(lapply(expo_INMA, function(x) log(x+0.0001)))
rownames(expoINMA.log) <- rownames(expo_INMA)

# =============================================================================
# EXAMPLE USAGE OF Fit() FUNCTION: Table of fits of first 100 probes:
# =============================================================================

source("Fit.R")

# Fit() attributes:
# _________________
# N = position (in exprs(expset))[N,] of probes which are going to be analysed
# select = global signifiance level for the MFP algorithm to be adjuted due to cor(expo)
# select_adj = manual adjusted select 
# cores = number of cores to be used (if > 1, parallel processes will be executed)

EXAMPLE  <- Fit(expset = expset_INMA, expo = expo_INMA, N = sample(1:nrow(exprs(expset_INMA)), 100), select = 0.05, cores = 20)

EXAMPLE  <- Fit(expset = expset_INMA, expo = expoINMA.log, N = sample(1:nrow(exprs(expset_INMA)), 100), df = 2, select = 0.05, cores = 20)


print(EXAMPLE[[1]]) # Table ordered by LRT sig.

print(EXAMPLE[[2]][1:10]) # List of associated expositions to the first 10 probes


# =============================================================================
# EXAMPLE USAGE OF SingleFit() FUNCTION: Fit 1 probe
# =============================================================================

source("SingleFit.R")

# SingleFit() attributes:
# _________________
# probe = name of position of probe (in exprs(expset)[probe,]) which is going to be analysed
# select = global signifiance level for the MFP algorithm to be adjuted due to cor(expo)
# select_adj = manual adjusted select (select attribute will be ignored)

example <- SingleFit(expset = expset_INMA, expo = expo_INMA, probe = "TC09001306.hg.1", df = 2, select = 0.05) # , select_adj = 0.05) 
# sig. probes = 3, TC01000083.hg.1, TC04000946.hg.1, TC12001765.hg.1, TC08001774.hg.1 
example <- SingleFit(expset = expset_INMA, expo = expoINMA.log, probe = "TC03002103.hg.1", df = 2, select = 0.05)



# summary.mfp.SingleFit() attributes:
# _________________
# mod = an mfp.SingleFit object
summary(example)


# plot.mfp.SingleFit() attributes:
# _________________
# mod = an mfp.SingleFit object
# realpoints = logical indicator if real points will be added to the plot of the marginal effects
# seed = seed to be fixed for the colors of the plot
pdf(paste0("plots/plot_", example$probe_name, ".pdf"), onefile = FALSE)
plot(example, realpoints = TRUE, seed = 2)
dev.off()




ff <- function(x, ...) {
 y <- names(x)
 example <- SingleFit(expset = expset_INMA, expo = expoINMA.log, 
                      probe = y, select = 0.05)
 plot(example, realpoints = TRUE, seed = 2, ...)   
 points(expoINMA.log[,grep(x, names(expoINMA.log))], exprs(expset_INMA)[y,], col="gray60", cex=0.7)
 lines(lowess(expoINMA.log[,grep(x, names(expoINMA.log))], exprs(expset_INMA)[y,]), col="lightblue")
 abline(reg=lm(expoINMA.log[,grep(x, names(expoINMA.log))] ~ exprs(expset_INMA)[y,]), col="red")
 dev.off()
}

ff(EXAMPLE[[2]][3], xlim=c(0,7), ylim=c(0,2))





