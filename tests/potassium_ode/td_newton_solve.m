function [W, Niter, Res] = td_newton_solve( G, J, W0, MaxIter, Tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A general purpose newton solver for vector functions 
% Target:
%   To find W in R^N, such that: G(W) = 0.
%   G : R^N -> R^N
% Usage:
%   [W, Niter, Res] = 
%     td_newton_solve( G, J, W0, MaxIter, Tol )
% Input:
%     G        : A function handle takes N-vec and 
%                returns N-vec
%     J        : A function handle takes N-vec and
%                returns an NxN matrix.
%     W0       : A N-vec as initial guess
%     MaxIter  : Maximum number of iterations
%     Tol      : Absolute tolerance
% Output:
%     W        : N-vec solution to the nonlinear eqn.
%     Niter    : Actual number of iterations
%     Res      : Residual of the final G(W)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  Niter = 0; % Iteration counter
  W     = W0;
  Gval  = G(W);
  Res   = norm(Gval);

  while (Niter < MaxIter && Res > Tol)
    Niter = Niter + 1;
    W     = W - inv(J(W))*Gval;
    Gval  = G(W);
    Res   = norm(Gval);
  end
  
end
