function u_next = solve_potassium_one_step_dirk2( opts, u_this, frhs, fjac )
  % Solving 
  %   du/dt = f(u) : fhs = f, fjac = j = df/du
  % using backward Euler:
  %   u_next - u_this = dt * f( u_next )
  % or equivalently F(u) = 0, u = u_next:
  %   F(u) = u - u_this - dt * f( u )
  % The Jacobian of F is:
  %   J(u) = I - dt * j( u )

  % L-stable DIRK2 (2nd order, Pareschi and Russo)
  a11 = 1.0-sqrt(2.0)/2.0; 
  a21 = 1.0-2.0*a11;
  a22 = a11;
  c1  = a11;
  c2  = a21 + a22;
  b1  = 1.0/2.0;
  b2  = 1.0/2.0;

  % Stage 1:
  % u_stg1 = u_this + a11 * dt * f( u_stg1 )
  a11dt = a11 * opts('dt');
  nt_F = @(u)u-u_this-a11dt*frhs(u);
  nt_J = @(u)eye(length(u))-a11dt*fjac(u);
  [u_stg1, iter, res] = td_newton_solve(nt_F, nt_J, u_this, opts('ntIter'), opts('ntTol'));
  if ( opts('show') ~= 0 )
    disp( sprintf('    DIRK2 - Stage 1 - Newton Report: %d iterations, residual = %e', iter, res) );
  end

  % Stage 2:
  % u_stg2 = u_this + a21 * dt * f( u_stg1 )
  %                 + a22 * dt * f( u_stg2 )
  f1    = frhs( u_stg1 );
  a21dt = a21 * opts('dt');
  f_fix = u_this + a21dt * f1;
  a22dt = a22 * opts('dt');
  nt_F = @(u)u-f_fix-a22dt*frhs(u);
  nt_J = @(u)eye(length(u))-a22dt*fjac(u);
  [u_stg2, iter, res] = td_newton_solve(nt_F, nt_J, u_this, opts('ntIter'), opts('ntTol'));
  if ( opts('show') ~= 0 )
    disp( sprintf('    DIRK2 - Stage 2 - Newton Report: %d iterations, residual = %e', iter, res) );
  end

  % Final solution
  u_next = u_this + b1 * opts('dt') * f1 + b2 * opts('dt') * frhs(u_stg2);
end
