function rhs = potassium_ode( u, pars )
  % u is an n x 1 vector; one should have n = 4
  rhs = zeros(size(u));
  B1 = u(1); B2 = u(2); P = u(3); V = u(4);
  H = 1.0/(1.0+exp(-2.0*pars('k')*(P-pars('Ps'))));
  rhs(1) = - pars('theta') * (1-H) * B1 * V;
  rhs(2) = pars('theta') * (1-H) * B1 * V - pars('delta') * B2;
  rhs(3) = pars('a') * pars('delta') *  B2 - pars('tau') * P;
  rhs(4) = pars('b') * pars('delta') * (1-H) * B2 - pars('gamma') * V - ...
           pars('theta') * (1-H) * B1 * V / pars('V0');
end
