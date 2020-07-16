prob_id       = -1
framework     = 3

lambda        = 0.0
mu            = 0.0
delta_n       = 0.0

r_0           = 1.0
Nx            = 50
cfl           = 0.8
Tinit         = 0.0
Tterm         = 2.0
Rterm         = 100.0
Rdead         = 0.01
Tfreq         = 0.5
//frequency     = 1
max_iter      = 100000
dt_init       = 1.0
dt_growth     = 1.2

lsysType      = DIR
lsysPrint     = 0

outputDataFile  = sol.dat
radiusDataFile  = rad.dat

init_u_par    = 0.5
n_bc          = 0.5
ic            = 0
time_method   = 0
space_cell    = 0
space_chem    = 0
dr_option     = 2
r_option      = 1

assemble          = 1
newton_iter       = 3
Jacobian_limiter  = 0
Jacobian_velocity = 1
chem_r_ppm        = 0
infiltration      = 4
