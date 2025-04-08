#author: Jamie Sanderlin
#date: 11/14/2024

#purpose: this script uses wrangled SW bird CASC data for spOoccupancy multi-species occupancy models
#         and runs the model

###################
# load functions  #
###################
source(here::here("Code","02-model","02-model.load.packages.R"))
#source(here::here("Code","02-model","02-model.inits.R"))

#######################################
# load data from data wrangling step  #
#######################################

wrangled.data <- readRDS(file=here::here("Output","firebird_modeldata_spOcc.rds"))

##############################################
# provide initial parameter values for model #
##############################################

inits.list <- list(alpha.comm = 0, 
                   beta.comm = 0, 
                   beta = 0, 
                   alpha = 0,
                   tau.sq.beta = 1, 
                   tau.sq.alpha = 1, 
                   z = apply(wrangled.data$Y, c(1, 2), max, na.rm = TRUE))

prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                   alpha.comm.normal = list(mean = 0, var = 2.72), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

#____________arguments___________________________
n.chains <- 3
n.samples <- 3000
n.thin <- 1
n.burn <- 2000

#parameters.to.save=c('omega','rho.za','zeta','beta','beta.transect','beta.year',
#                     'beta.veg','zeta.obs')

occ.formula = ~ simpson.nfire*fire.ind.transect + sev.SD*fire.SD.ind.transect+
  resid.tsf.sev*fire.ind.pnt+resid.tsf2.sev2*fire.ind.pnt+PC1+PC2+fcv+TRI+hli+stream+elev+lat+
  (1|transectID)+(1|yearID)+(1|vegtype)

gc()

out <- msPGOcc(occ.formula = occ.formula, 
               det.formula = ~ DOY + tday + countlength + do.ind, 
               data = wrangled.data, 
               inits = inits.list, 
               n.samples = n.samples, 
               priors = prior.list, 
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = n.burn, 
               n.thin = n.thin, 
               n.chains = n.chains)

summary(out, level = 'community')
summary(out, level = 'species')
