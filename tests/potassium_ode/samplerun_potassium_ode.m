setup_potassium_prob;

% B1, B2, P, V
u0 = [1.0; 0.0; 0.0; 1.0];

[time, u] = solve_potassium_ode( u0, pars, opts );
