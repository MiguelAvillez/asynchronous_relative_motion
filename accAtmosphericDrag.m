%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion Approximations 
%   Under J2 and Atmospheric Drag", TODO
%
% Summary:
%   Computes the acceleration wrt RTN frame acting on an individual satellite due
%     to the atmospheric drag perturbation.
%   Assumptions: cannonball drag model, static atmosphere
%
% Inputs:
%   tt: argument of latitude
%   state: absolute state of the spacecraft: [bb; x; y; p; oo; t]
%   mu: gravitational parameter
%   R: Radius of the central planet
%   j2: J2 coefficient of the gravity model
%   k: inverse of the ballistic parameter
%   rho: atmospheric density
%
% Outputs:
%   acc: acceleration in RTN frame: [ar; af; ah]
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: August 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acc = accAtmosphericDrag(tt, state, mu, R, j2, k, rho)

  % Extract state elements
  bb = state(1);
  x = state(2);
  y = state(3);

  % Remove normalization from eccentricities
  ex = j2 * x;
  ey = j2 * y;

  % Compute RTN acceleration
  ar = (1/2).*bb.^2.*k.*mu.*R.^(-1).*rho.*(ey.*cos(tt)+(-1) ...
    .*ex.*sin(tt)).*(1+ex.^2+ey.^2+2.*ex.*cos(tt)+2.*ey.*sin(tt)).^( ...
    1/2);
  af = (-1/2).*bb.^2.*k.*mu.*R.^(-1).*rho.*(1+ex.*cos(tt)+ ...
    ey.*sin(tt)).*(1+ex.^2+ey.^2+2.*ex.*cos(tt)+2.*ey.*sin(tt)).^(1/2);
  ah = 0;

  % Concatenate acceleration components
  acc = [ar; af; ah];

end