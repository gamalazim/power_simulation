##
## G. Abdel-Azim
## Sample size simulation for a categorical variable with 2 factors 
## November 2018
##

# Function: getThresholds
# Return two thresholds to divide a continuous effect into 3 categories 
getThresholds <- function(cont.eff, group.diff, min.cat = 100) {
  loss <- 100
  thr <- numeric(2)
  thr[1] <- summary(cont.eff)['1st Qu.']
  thr[2] <- summary(cont.eff)['3rd Qu.']
  eps <- mean(c(thr[1], thr[2]))/20
  t.change <- 2; if(runif(1) <= 0.5) t.change <- 1; 
  
  while(loss > 0.1) {
    t.change.old <- t.change; thr.old <- thr
    t.change <- 2; if(runif(1) <= 0.5) t.change <- 1;
    cat.len <- c(length(cont.eff[cont.eff<thr[1]]), 
                 length(cont.eff[cont.eff >= thr[1] & cont.eff < thr[2] ]),
                 length(cont.eff[cont.eff>=thr[2]]))
    if( any(cat.len) < min.cat) { thr <- thr.old; break }
    
    mean.diff <- mean(cont.eff[cont.eff >= thr[2]]) - mean(cont.eff[cont.eff < thr[1]])
    
    if(mean.diff <= 0) { thr <- thr.old; break }
    group.diff.diff <- abs(mean.diff - group.diff)
    number.diff <- abs(length(cont.eff[cont.eff >= thr[2]]) - length(cont.eff[cont.eff < thr[1]]))
    
    loss <- group.diff.diff * number.diff
    
    if(mean.diff < group.diff) ifelse(t.change==1, thr[1] <- thr[1]-eps, thr[2] <- thr[2]+eps)
    if(mean.diff > group.diff) ifelse(t.change==1, thr[1] <- thr[1]+eps, thr[2] <- thr[2]-eps)
  }
  return(thr)
}

#
# Main simulation
#

nrep <- 100  # Assign number of replicates
REP <- matrix(0, nrow=nrep, ncol=5) ## Results matrix
ii <- 0 # index of repeat, don't change! 
repeat {
  ii <- ii+1; cat("Simulation Replicate: ", ii, "\n")
  Ns <- 12 ; sires <- data.frame(sire=1:Ns, effects = runif(Ns, .3, .35), no_batches = sample(2:8, size=Ns, replace=T))
  Nb <- sum(sires$no_batches); 
  batches <- data.frame(batch=1:Nb, effects = runif(Nb, 0, .03), no_obs = sample(500:700, size=Nb, replace=T)) #100-200
  N <- sum(batches$no_obs)
  batch_cumsum <- c(0, cumsum(sires$no_batches))
  obs_sum <- numeric(Ns)
  for(i in 1:Ns) obs_sum[i] <- sum(batches[(batch_cumsum[i]+1) : batch_cumsum[i+1], "no_obs"])
  
  d <- data.frame(
    batch = factor(rep(batches$batch, batches$no_obs)),
    batch_effect = rep(batches$effects, batches$no_obs),
    sire = factor(rep(sires$sire, obs_sum)),
    sire_effect = rep(sires$effects, obs_sum)
  )
  
  #
  # Create 2 extreme categories with a difference of 'group.diff' between category averages 
  #
  group.diff <- 0.02
  thr <- getThresholds(d$batch_effect, group.diff)
  
  d$vgroup <- 0
  d[d$batch_effect > thr[1], "vgroup"] <- 1
  d[d$batch_effect >= thr[2], "vgroup"] <- 2
  d$vgroup <- factor(d$vgroup)
  #cat("group sizes: ", table(d$vgroup), " -- ")
  #diff(tapply(d$batch_effect, d$vgroup, mean), lag = 2)
  
  if(length(table(d$vgroup)) != 3) { ii <- ii - 1; next }
  
  # add up batch_effect + sire_effect + residuals to form the linear response variable
  d$y <- d$batch_effect + d$sire_effect + rnorm(N, 0, .015)
  d$batch <- factor(d$batch)
  d$sire <- factor(d$sire)
  d$z <- numeric(N)
  for(i in 1:N) d$z[i] <- sample(0:1, prob=c((1-d$y[i]), d$y[i]), size = 1)
  # run a glm with bull + group - binomial family with default link 
  d.glm <- glm(z ~ vgroup + sire, d, family = binomial)
  REP[ii,] <- c(summary(d.glm)$coeff[3,c(1,4)], N, mean(table(d$vgroup)[c(1,3)]), 
                diff(tapply(d$batch_effect, d$vgroup, mean), lag = 2))
  colnames(REP) <- c("Est", "PVal", "N", "AvgGrpSz", "GrpDiff")
  
  if(ii == nrep) break
} ## End repeat


##
## Simulation summary
##
cat("\n ================================================\n", 
    "Results averaged over all simulation replicates:\n",
    "================================================\n",
    "POWER (alpha = 5%):             ", nrow(REP[REP[,"PVal"] <= .05,])/nrep, "\n",
    "Liner scale difference:         ", mean(REP[,"Est"]), "\n",
    "Total data size:                ", mean(REP[,"N"]), "\n",
    "Avg no. of 2 extreme groups:    ", mean(REP[,"AvgGrpSz"], na.rm = T), "\n",
    "True extreme group difference:  ", mean(REP[,"GrpDiff"], na.rm = T), "\n")

## Done









