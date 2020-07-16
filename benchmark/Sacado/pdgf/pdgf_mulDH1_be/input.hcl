prob_id       = 0
framework     = 3

lambda        = 0.48
alpha         = 0.06
rho           = 0.9
gamma         = 100.0
delta         = 0.45
mu            = 0.33
D             = 6.048
theta         = 1.0e+6
beta          = 1.0e+5
m             = 1.5e+5

r_0           = 0.2
Nx            = 050
cfl           = 8.0
Tinit         = 0.0
Tterm         = 100.0
//Rterm         = 0.5
Rterm         = 5.5
Rdead         = 0.01
Tfreq         = 0.02
//frequency     = 1
max_iter      = 10000000
dt_init       = 1.0e-6
dt_growth     = 1.2

outputDataFile  = sol.dat
radiusDataFile  = rad.dat

lsysType      = DIR
lsysPrint     = 0

A_Multiplier  = 2
n_bc          = 5.0e-3
ic            = 0
time_method   = 0
space_cell    = 0
space_chem    = 0
dr_option     = 2
r_option      = 1

assemble          = 1
newton_iter       = 2
Jacobian_limiter  = 0
Jacobian_velocity = 1
chem_r_ppm        = 0
