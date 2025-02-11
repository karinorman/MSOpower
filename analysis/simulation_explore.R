library(spOccupancy)

### Fixed Parameters ###
## effect_size = .25 decline
## study_duration = 10
## single deployment

### Parameters ###
psi1 = 0.43 #initial ocucpancy, from Wood 2019
p = c(0.4, 0.8) #detection probability, from Wood 2019
epsilon1 = .2 #extinctino probability, Wood 2019, based on biology?

# derive from psi1 = gamma/(gamma + epsilon); occupancy at equilibrium
gamma1 = (psi1 * epsilon1)/(1-psi1)


### Design scenarios ###
n_sites
n_recorders


### Simulate decline by multiplicative decline in survival

test_data <- simIntOcc(J.x = 100, J.y = 100,
          n.rep = list(rep(2, 100*100)), beta = c(psi1, 0.2), alpha = list(c(p[1], 0.2)))

# add trend
trend_data <- simTIntOcc()
