
# load libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# reactor parameters
L = 10
vz = 1
D = 1e-4
D_uL = D/(vz*L)
tau = L/vz

# reaction parameters
# 1. A + B -> C
# 2. B + C -> D
#
k1 = 1 
k2 = 1



# reactor inlet conditions
cinlet = c(1.0, 1.0, 0.0, 0.0)

# other parameters
ngrid = 100
nc = 4 

# fluctuation parameters for cB_inlet;
incFluct = 0
tstart = 1
amp_pct = 0.5
amp_freq = 0.1 


parms = list(D_uL = D_uL, ngrid = 100, nc = 4, cinlet = cinlet,k1 = k1*tau, k2 = k2*tau,
             tstart = tstart, amp_pct = amp_pct, amp_freq = amp_freq, incFluct = 1
             ) 

cinit = rep(0,parms$ngrid*parms$nc)
tspan = seq(0, 3, by = 0.01)
sol =ode(y = cinit, times = tspan, func = ratefn, parms = parms, jactype = "bandint",
            bandup = 12, banddown = 12)

cf = sol[,c(1,(parms$ngrid-1)*parms$nc + c(2,3,4,5))]
cf = data.frame(cf)
names(cf) = c("time","cA","cB","cC","cD")
cf_g = gather(cf,var,value,-time)

ntime = nrow(sol)
cfL = sol[ntime,-1]
dz = 1/(ngrid + 1)
zL = rep(seq(dz,ngrid*dz, by = dz), each = nc)
compname = rep(c("A","B","C","D"),ngrid)

cfL_df = data.frame(compname = compname, zL = zL, cfL = cfL)

ggplot(cf_g) + geom_line(aes(x = time, y = value, color = var)) + theme_bw()

ggplot(cfL_df) + geom_line(aes(x = zL, y = cfL, color = compname)) + theme_bw()
