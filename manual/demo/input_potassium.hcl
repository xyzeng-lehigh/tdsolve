// Example input file for the potassium model (see manual/doc/potassium.tex)
//   TumorDynRun -i input_potassium.hcl
// Lines that start with // or # are comments.

// ---- Model and geometry -----------------------------------------------
// potassium model (3), solved by the generic framework (3)
prob_id   = 3
framework = 3
// 1: slab, 2: radial symmetry in the plane, 3: spherical symmetry
dimension = 2
// radius R of the (fixed) domain
r_0       = 100.0
// number of cells per unit volume
theta     = 2.96e+5

// ---- Parameters of the model ------------------------------------------
potassium
  global
    // ratio rho between the intracellular and the extracellular space
    ie_ratio       = 10.0
    // 1: release of P and V with a time delay, 0: without
    enable_delay   = 1.0
  end

  B1
    // infection rate beta
    grow_neg_virus = 0.5e-9
    // fraction c_1 of the infection that is blocked by potassium
    grow_block     = 0.75
  end

  B2
    // death rate delta of the infected cells
    grow_neg       = 0.08
    // fraction c_2 of the death rate that is blocked by potassium
    grow_block     = 0.5
  end

  P
    // decay rate tau and diffusion coefficient D_P
    grow_neg       = 0.0
    diffusion      = 1.0e-1
    // steepness k and threshold P_s of the blocking function
    slope          = 4.0
    threshold      = 25.0
    // potassium released by a dying cell (a)
    burst          = 5.0625e-4
  end

  V
    // decay rate gamma and diffusion coefficient D_V
    grow_neg       = 0.0
    diffusion      = 1.0e-2
    // number of viruses N_V placed at the origin at the initial time
    init_num       = 1.0e+5
    // viruses released by a dying cell (b), and the blocked fraction c_V
    burst          = 3.0e+2
    burst_block    = 0.5
  end
end

// ---- Mesh and time stepping -------------------------------------------
// number of mesh intervals in [0, R]
Nx             = 100
// initial and final time
Tinit          = 0.0
Tterm          = 400.0
// first time step, largest and smallest time step, growth per step
dt_init        = 0.2
dt_max         = 0.2
dt_min         = 0.2
dt_growth      = 1.2
// largest number of time steps
max_iter       = 10000000

// ---- Numerical method -------------------------------------------------
// 1: implicit (backward Euler with Newton's method), 0: explicit for B1-B3
assemble       = 1
// keep 0 for this model when enable_delay = 1
time_method    = 0
// largest number of Newton iterations per time step, and their tolerance
newton_iter    = 100
newton_tol     = 1.0e-8
// direct solver for the linear systems, without printing their residuals
lsysType       = DIR
lsysPrint      = 0

// ---- Output -----------------------------------------------------------
// time between two outputs of the solution
Tfreq           = 100.0
outputDataFile  = sol.dat
radiusDataFile  = rad.dat
