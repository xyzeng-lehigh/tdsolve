function jac = potassium_ode_jac( u, pars )
  % u is an n x 1 vector; one should have n = 4
  jac = zeros(length(u));
  B1 = u(1); B2 = u(2); P = u(3); V = u(4);
  H = 1.0/(1.0+exp(-2.0*pars('k')*(P-pars('Ps'))));
  dHdP = 2.0*pars('k')*exp(-2.0*pars('k')*(P-pars('Ps')))*H*H;
  jac(1,1) = - pars('theta') * (1-H) * V;
  jac(1,3) = - pars('theta') * (-dHdP) * B1 * V;
  jac(1,4) = - pars('theta') * (1-H) * B1;
  jac(2,1) = pars('theta') * (1-H) * V;
  jac(2,2) = - pars('delta');
  jac(2,3) = pars('theta') * (-dHdP) * B1 * V;
  jac(2,4) = pars('theta') * (1-H) * B1;
  jac(3,2) = pars('a') * pars('delta'); 
  jac(3,3) = - pars('tau');
  jac(4,1) = - pars('theta') * (1-H) * V / pars('V0');
  jac(4,2) = pars('b') * pars('delta') * (1-H);
  jac(4,3) = pars('b') * pars('delta') * (-dHdP) * B2 - ...
             pars('theta') * (-dHdP) * B1 * V / pars('V0');
  jac(4,4) = - pars('theta') * (1-H) * B1 / pars('V0') - pars('gamma');
end
