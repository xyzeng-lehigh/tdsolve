prob_id   = 3
framework = 3
theta     = 2.96e+5
dimension = 2

potassium
  global
    ie_ratio       = 10.0
    enable_delay   = 1.0
  end

  B1
    grow_neg_virus = 0.5e-9
    grow_block     = 0.75
  end

  B2
    grow_neg       = 0.08
    grow_block     = 0.5
  end
  
  P
    grow_neg       = 0.0
    diffusion      = 1.0e-1
    slope          = 4.0
    threshold      = 25.0
    burst          = 5.0625e-4
  end

  V
    diffusion      = 1.0e-2
    init_num       = 1.0e+5
    burst          = 3.0e+2
    burst_block    = 0.5
  end
end

r_0            = 100.0
cfl            = 1.0
Nx             = 10
Tinit          = 0.0
Tterm          = 100.0
Rterm          = 100.0
Rdead          = 0.1
Tfreq          = 10.0
max_iter       = 10000000
dt_init        = 0.2
dt_max         = 0.2
dt_min         = 0.2
dt_growth      = 1.2
auto_reduce_dt = 1

outputDataFile  = sol.dat
radiusDataFile  = rad.dat

lsysType      = DIR
lsysPrint     = 0

ic            = 0
time_method   = 0
space_cell    = 0
space_chem    = 0

assemble          = 1
newton_iter       = 100
newton_tol        = 1.0e-8
