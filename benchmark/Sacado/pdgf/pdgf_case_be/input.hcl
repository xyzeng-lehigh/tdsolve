prob_id       = 0
framework     = 3

lambda        = 1.0
alpha         = 1.0
rho           = 0.0
gamma         = 0.0
delta         = 0.0
mu            = 0.0
D             = 1.0
theta         = 1.0
beta          = 1.0
m             = 30.0

r_0           = 1.0
Nx            = 050
cfl           = 8.0
Tinit         = 0.0
Tterm         = 1.0
Rterm         = 100.0
Rdead         = 0.01
Tfreq         = 0.2
//frequency     = 1
max_iter      = 100000
dt_init       = 1.0
dt_growth     = 1.2
auto_reduce_dt = 0

outputDataFile  = sol.dat
radiusDataFile  = rad.dat

lsysType      = DIR
lsysPrint     = 0

A_Multiplier  = 2
n_bc          = 0.5
ic            = 2
time_method   = 0
space_cell    = 0
space_chem    = 0
cell_rflux    = 1
chem_rflux    = 1
dr_option     = 2
r_option      = 1

assemble          = 1
newton_iter       = 3
Jacobian_limiter  = 0
Jacobian_velocity = 1
