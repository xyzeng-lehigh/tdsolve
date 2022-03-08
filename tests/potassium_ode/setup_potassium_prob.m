keys = {'N', 'theta', 'delta', 'a', 'b', 'gamma', 'tau', 'k', 'Ps', 'V0'};
vals = {0.0, 0.0    , 0.0    , 0.0, 0.0, 0.0    , 0.0  , 0.0, 0.0 , 0.0 };
pars = containers.Map( keys, vals ); % Coefficients for the non-dimensional ODE
refs = containers.Map( keys, vals ); % Reference values

refs('N')     = 3.0e+4;
refs('theta') = 1.5e-8;
refs('delta') = 0.06;
refs('a')     = 4.5e+1;
refs('gamma') = 0.004;
refs('tau')   = 0.004;
refs('k')     = 1.0e-4;
refs('Ps')    = 1.0e+6;
refs('V0')    = 1.0e+4;
refs('b')     = 1.0; 

pars('N')     = 1.0; % Not to be used at all
pars('theta') = refs('theta') * refs('V0'); 
pars('delta') = refs('delta');
pars('a')     = refs('a') * refs('N') / refs('Ps');
pars('gamma') = refs('gamma');
pars('tau')   = refs('tau');
pars('k')     = refs('k') * refs('Ps');
pars('Ps')    = 1.0;
pars('V0')    = refs('V0') / refs('N');
% A major parameter to vary
pars('b')     = refs('b') * refs('N') / refs('V0');


s_keys = {'dt', 'T', 'ntIter', 'ntTol', 'show', 'time'};
s_vals = {0.0,  0.0, 100     , 1e-10  , 1     , 'be'  };
opts   = containers.Map( s_keys, s_vals );

opts('dt')     = 1.0;
opts('T')      = 4800.0;
opts('ntIter') = 100;
opts('ntTol')  = 1e-8;
opts('show')   = 1;
opts('time')   = 'dirk2'; % 'be' or 'dirk2'

% Compute rho using dimensionless parameters
frho = @(theta,b,V0,gamma,k)theta*b/(theta/V0+gamma+(theta/V0+2.0*gamma)*exp(-2*k)+gamma*exp(-4*k));
fb   = @(theta,rho,V0,gamma,k)rho*(theta/V0+gamma+(theta/V0+2.0*gamma)*exp(-2*k)+gamma*exp(-4*k))/theta;
