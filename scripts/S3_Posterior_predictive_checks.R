library(spOccupancy)

load("models/22222/sfMsPGOcc/model_sfMsPGOcc_1981-2010_nthin50_nbatch10000_nchain4_nburn150000.RData")
load("models/22222/spPGOcc/model_spPGOcc_1981-2010_nthin500_nbatch1e+05_nchain4_nburn1500000.RData")


# [sfMsPGOcc] Posterior Predictive Checks ---------------------------------------------
# Freeman-tukey 
ppc.sfMsPGOcc.ft1 <- ppcOcc(out.sfMsPGOcc, fit.stat = 'freeman-tukey', group = 1) # group by sites
ppc.sfMsPGOcc.ft2 <- ppcOcc(out.sfMsPGOcc, fit.stat = 'freeman-tukey', group = 2) # group by replicates

summary(ppc.sfMsPGOcc.ft1)
summary(ppc.sfMsPGOcc.ft2)

# Chi-squared
ppc.sfMsPGOcc.cq1 <- ppcOcc(out.sfMsPGOcc, fit.stat = 'chi-squared', group = 1) # group by sites
ppc.sfMsPGOcc.cq2 <- ppcOcc(out.sfMsPGOcc, fit.stat = 'chi-squared', group = 2) # group by replicates

summary(ppc.sfMsPGOcc.cq1)
summary(ppc.sfMsPGOcc.cq2)



# [spPGOcc] Posterior Predictive Checks -----------------------------------
# Freeman-tukey 
ppc.spPGOcc.ft1 <- ppcOcc(out.spPGOcc, fit.stat = 'freeman-tukey', group = 1) # group by sites
ppc.spPGOcc.ft2 <- ppcOcc(out.spPGOcc, fit.stat = 'freeman-tukey', group = 2) # group by replicates

summary(ppc.spPGOcc.ft1)
summary(ppc.spPGOcc.ft2)

# Chi-squared
ppc.spPGOcc.cq1 <- ppcOcc(out.spPGOcc, fit.stat = 'chi-squared', group = 1) # group by sites
ppc.spPGOcc.cq2 <- ppcOcc(out.spPGOcc, fit.stat = 'chi-squared', group = 2) # group by replicates

summary(ppc.spPGOcc.cq1)
summary(ppc.spPGOcc.cq2)


# WAIC --------------------------------------------------------------------
waic.sfMsPGOcc = waicOcc(out.sfMsPGOcc)
waic.spPGOcc = waicOcc(out.spPGOcc)

load("output/ppc_results.RData")



if(!dir.exists("output")){
  dir.create("output")
}

# Save the posterior predictive check results
save(
  # Multi-species model PPC results
  ppc.sfMsPGOcc.ft1, ppc.sfMsPGOcc.ft2,
  ppc.sfMsPGOcc.cq1, ppc.sfMsPGOcc.cq2,
  
  # Single-species model PPC results
  ppc.spPGOcc.ft1, ppc.spPGOcc.ft2,
  ppc.spPGOcc.cq1, ppc.spPGOcc.cq2,  
  
  # WAIC results
  waic.sfMsPGOcc, waic.spPGOcc,
  
  file = "output/ppc_results.RData"
)

