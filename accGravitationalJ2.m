%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion Approximations 
%   Under J2 and Atmospheric Drag", TODO
%
% Summary:
%   Computes the acceleration wrt RTN frame acting on an individual satellite due
%   to the J2 gravitational perturbation.
%
% Inputs:
%   tt: argument of latitude
%   state: absolute state of the spacecraft: [bb; x; y; p; oo; t]
%   mu: gravitational parameter
%   R: Radius of the central planet
%   j2: J2 coefficient of the gravity model
%
% Outputs:
%   acc: acceleration in RTN frame: [ar; af; ah]
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: August 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acc = accGravitationalJ2(tt, state, mu, R, j2)

  % Extract state elements
  bb = state(1);
  x = state(2);
  y = state(3);
  p = state(4);

  % Remove normalization from eccentricities
  ex = j2 * x;
  ey = j2 * y;

  % Compute RTN acceleration
  ar = (-3/2).*bb.^8.*j2.*mu.*R.^(-2).*(1+ex.*cos(tt)+ey.*sin(tt)).^4.*( ...
    1+3.*((-1)+bb.^2.*p.^2).*sin(tt).^2);
  af = 3.*bb.^8.*j2.*mu.*((-1)+bb.^2.*p.^2).*R.^(-2).*cos(tt).*sin(tt).*( ...
    1+ex.*cos(tt)+ey.*sin(tt)).^4;
  ah = (-3).*bb.^9.*j2.*mu.*p.*(1+(-1).*bb.^2.*p.^2).^(1/2).*R.^(-2).* ...
    sin(tt).*(1+ex.*cos(tt)+ey.*sin(tt)).^4;

  % Concatenate acceleration components
  acc = [ar; af; ah];

end