prob_id   = 2
framework = 3
theta     = 1.0e+6

oncolyticOld
  B
    init = 0.999995
    mm_ceil_grow_pos    = 0.0120
    mm_half             = 317
    coef_vi_neg         = 0.89e-9
    coef_h_neg          = 3.42e-8
    mm_ceil_coef_mb_neg = 0.0
  end

  H
    init = 1e-6
    mm_ceil_grow_pos    = 0.01
    mm_half             = 317
    mm_ceil_coef_mb_neg = 0.0
    mm_ceil_coef_mv_neg = 0.0
    grow_neg            = 2.35e-3
    prod_vi_b_pos       = 0.89e-9
    coef_b_pos          = 3.42e-8
  end

  N
    init = 1e-6
    mm_half               = 317
    mm_ceil_prod_mb_b_pos = 0.0
    mm_ceil_prob_mb_h_pos = 0.0
    mm_ceil_prod_mv_h_pos = 0.0
    prod_h_pos            = 2.35e-3
  end

  G
    init = 1e-6
    coef_b_pos = 1.37e-8
    grow_neg   = 0.03125
  end

  Mv
    init = 1e-6
    mm_ceil_grow_pos = 0.02
    mm_half          = 317
    grow_neg         = 0.2
    coef_h_pos       = 5.6e-7
  end

  Mb
    init = 1e-6
    mm_ceil_grow_pos = 0.02
    mm_half          = 317
    grow_neg         = 0.2
    coef_b_pos       = 5.6e-7
    coef_h_pos       = 5.6e-7
  end

  A
    init = 0.0
    mm_half             = 317
    coef_mm_b_neg       = 5.255e-6
    mm_ceil_b           = 0.0120
    coef_mm_h_neg       = 7.2e-6
    mm_ceil_h           = 0.01
    coef_mm_mv_neg      = 1.5e-6
    mm_ceil_mv          = 0.02
    coef_mm_mb_neg      = 1.5e-6
    mm_ceil_mb          = 0.02
    prod_g_neg          = 0.00375
    prod_constant       = 30
    diffusion           = 2.42
  end

  Vi
    init = 0.0
    prod_growth_h_pos = 2.5
    growth_h          = 2.35e-3
    grow_neg          = 1.85e-3
    diffusion         = 2.25e-2
  end
end

r_0           = 0.2
cfl           = 8.0
Nx            = 40
Tinit         = 0.0
Tterm         = 800.0
Rterm         = 0.22
Rdead         = 0.1
Tfreq         = 4.0
max_iter      = 10000000
dt_init       = 1.0e-6
dt_growth     = 1.2

outputDataFile  = sol.dat
radiusDataFile  = rad.dat

lsysType      = DIR
lsysPrint     = 0

ic            = 0
time_method   = 0
space_cell    = 0
space_chem    = 0
cell_rflux    = 1
chem_rflux    = 1
dr_option     = 2
r_option      = 1

posit_fix         = 1
assemble          = 1
newton_iter       = 3
Jacobian_limiter  = 0
Jacobian_velocity = 1
