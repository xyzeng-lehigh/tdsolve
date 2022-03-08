function [time, u] = solve_potassium_ode( u0, pars, opts )
  cur    = 0.0;
  time   = [cur];
  u      = u0;
  u_this = u0;
  iter   = 0;
  while cur < opts('T') 
    frhs = @(var)potassium_ode(var, pars);
    fjac = @(var)potassium_ode_jac(var, pars);
    if strcmp( opts('time'), 'be' )
      u_next = solve_potassium_one_step_be( opts, u_this, frhs, fjac );
    elseif strcmp( opts('time'), 'dirk2' )
      u_next = solve_potassium_one_step_dirk2( opts, u_this, frhs, fjac );
    end
    time   = [time cur];
    u      = [u u_next];
    u_this = u_next;
    cur = cur + opts('dt');
    iter = iter + 1;
    disp( sprintf( 'Iter %d, time %f', iter, cur ) );
  end
end
