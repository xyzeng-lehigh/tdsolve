function u_next = solve_potassium_one_step_be( opts, u_this, frhs, fjac )
  % Solving 
  %   du/dt = f(u) : fhs = f, fjac = j = df/du
  % using backward Euler:
  %   u_next - u_this = dt * f( u_next )
  % or equivalently F(u) = 0, u = u_next:
  %   F(u) = u - u_this - dt * f( u )
  % The Jacobian of F is:
  %   J(u) = I - dt * j( u )
  nt_F = @(u)u-u_this-opts('dt')*frhs(u);
  nt_J = @(u)eye(length(u))-opts('dt')*fjac(u);
  [u_next, iter, res] = td_newton_solve(nt_F, nt_J, u_this, opts('ntIter'), opts('ntTol'));
  if ( opts('show') ~= 0 )
    disp( sprintf('    BE - Newton Report: %d iterations, residual = %e', iter, res) );
  end
end
