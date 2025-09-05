%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion Approximations 
%   Under J2 and Atmospheric Drag", TODO
%
% Summary: Applies the solutions from section V.A and V.B of the paper
%   Uses the following analytical transformations:
%   - Osculating to mean relative state
%   - Mean to osculating relative state
%
% Note: the script computes the numerical mean relative state history via numerical quadrature. This computation
%       becomes a bit slow when sampling very large number of points. It is possible to deal with this by either
%       increasing the tolerance of the numerical quadrature or reducing the number of points.
%
% Dependencies:
%   oscToMeanRelativeState.m
%   meanToOscRelativeState.m
%   accAtmosphericDrag.m
%   accGravitationalJ2.m
%   stateArgLatDerivative.m
%
% Inputs:
%   Constants:
%       R: Radius of the central planet
%       mu: gravitational parameter
%       j2: J2 coefficient of the gravity model
%       we: Earth's angular velocity
%   Initial argument of latitude: tti
%   Initial osculating Keplerian orbital elements of the chief:
%       sma: semi-major axis
%       ex: x-eccentricity
%       ey: y-eccentricity
%       inc: inclination
%       raan: right ascension of ascending node
%       t: time
%   Initial osculating relative state of the chief in proposed orbital elements
%       dbb: relative beta
%       dx: relative X
%       dy: relative Y
%       dp: relative p
%       doo: relative right ascension of ascending node
%       dt: relative time
%   Drag parameters of the chief:
%       kC: ballistic coefficient
%       rhoC: atmospheric density
%   Drag parameters of the deputy:
%       kD: ballistic coefficient
%       rhoD: atmospheric density
%   Script variables:
%       numRevolutions: number of revolutions to compute
%       numPointsPerRev: number of points per computed revolution
%       expansionOrder: order of the power expansion -> allowed values are 2 and 3
%
% Outputs:
%   Plot of osculating->mean elements transformation
%   Plot of mean->osculating elements transformation
%   
%
% Authors: Miguel Avillez and David Arnas
% Created: August 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% Define parameters and select desired order of expansion

% Define constants
R = 6378.15e3; % central planet radius, m
mu = 3.986e14; % gravitational parameter, m^3/s^2
j2 = 1.08262668e-3; % J2 coefficient of the gravity model, -
we = 7.2921151467e-5; % Earth's angular velocity, rad/s
deg2rad = pi/180;
rad2deg = 1/deg2rad;

% Initial (i) argument of latitude (tt)
tti = 0 * deg2rad; % rad

% Define chief's orbit using Keplerian elements
% Sun-synchronous frozen orbit
smaC = 7077.722e3; % semi-major axis, m
exC = 4.5742e-4; % x-eccentricity (ex), - 
eyC = 0; % y-eccentricity (ey), -
incC = 98.186 * deg2rad; % inclination, rad
raanC = 0 * deg2rad; % right ascension of ascending node, rad
tC = 0; % time, s

% Define drag properties of chief
cdC = 2.2; % drag coefficient, -
sC = 8.5; % surface area, m2
mC = 1285; % mass, kg
kC = cdC * sC / mC; % inverse of ballistic coefficient
rhoC = 2e-13; % atmospheric density, kg/m3

% Define initial relative state using proposed orbital elements: [dbb, dx, dy, dp, doo, dt]
relativeStateInitial = [...
    -6.123e-06; % beta, -
    j2/2; % X, -,
    j2/5; % Y, -,
    -j2/10; % p, -
    pi; % right ascension of ascending node, rad
    -3e3]; % time, s

% Define drag properties of deputy
rhoD = rhoC; % atmospheric density, kg/m3
kD = 2 * kC; % inverse of ballistic coefficient

% Define auxiliary variable
dzeta = kD * rhoD - kC * rhoC;

% Select number of revolutions and number of points used for plotting
numRevolutions = 5;
if numRevolutions <= 10
    numPointsPerRev = 50;
else
    numPointsPerRev = 10;
end

% Order of the expansion. Allowed values: 2, 3
expansionOrder = 3;

%% Compute initial state of the chief and deputy with proposed orbital elements

betaC = sqrt(R/(smaC * (1-exC^2-eyC^2)));
xC = exC/j2;
yC = eyC/j2;
pC = cos(incC) / betaC;

% Initial state of the chief
stateInitialChief = [betaC; xC; yC; pC; raanC; tC];

% Initial state of the deputy
stateInitialDeputy = stateInitialChief + relativeStateInitial;

%% Numerically compute osculating relative state history

% Define the history of the independent variable (argument of latitude)
numPoints = numPointsPerRev * numRevolutions;
ttHistory = linspace(tti, numRevolutions * 2*pi, numPoints);

% Define numerical integration settings
tolerance = 1e-13;
odeOptions = odeset('RelTol', tolerance, 'AbsTol', tolerance);

% Define acceleration functions for chief (C) and deputy (D)
accFunctionC = @(tt,x) (accGravitationalJ2(tt, x, mu, R, j2) + accAtmosphericDrag(tt, x, mu, R, j2, kC, rhoC));
accFunctionD = @(tt,x) (accGravitationalJ2(tt, x, mu, R, j2) + accAtmosphericDrag(tt, x, mu, R, j2, kD, rhoD));

% Numerically propagate dynamics of chief and deputy
[~, propagatedStateHistC] = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionC),...
    ttHistory, stateInitialChief, odeOptions);
[~, propagatedStateHistD] = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionD),...
     ttHistory, stateInitialDeputy, odeOptions);

% Get osculating relative state history
propagatedOsculatingRelStateHist = propagatedStateHistD - propagatedStateHistC;

%% Compute mean state by applying the transformation osculating->mean to numerical propagation

% We take the numerical relative and absolute state histories: relOscState(tt) and oscStateC(tt)
% Then for each value of tt, we transform: [relOscState(tt), oscStateC(tt), tt] -> relMeanState(tt)

osc2meanAnalyticalRelStateHist = zeros(length(ttHistory),6);
for i = 1:length(ttHistory)

    % Apply analytical transformation from osculating to mean elements
    osc2meanAnalyticalRelStateHist(i,:) = oscToMeanRelativeState(...
        ttHistory(i), mu, j2, R, propagatedOsculatingRelStateHist(i,:), propagatedStateHistC(i,:), kD, rhoD, kC, rhoC, expansionOrder);
        
end

%% Numerically compute mean relative state history

% Numerically propagate full dynamics: backward for pi
[~, auxStateHistD] = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionD),...
    [ttHistory(1), ttHistory(1)-pi], stateInitialDeputy, odeOptions);
[~, auxStateHistC] = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionC),...
    [ttHistory(1), ttHistory(1)-pi], stateInitialChief, odeOptions);

% Propagate dynamics over range [ttHistory(1)-pi, ttHistory(end)+pi], which is range necessary to 
% evaluate the history of the mean relative state
stateHistDictD = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionD),...
    [ttHistory(1)-pi, ttHistory(end)+pi], auxStateHistD(end,:), odeOptions);
stateHistDictC = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionC),...
    [ttHistory(1)-pi, ttHistory(end)+pi], auxStateHistC(end,:), odeOptions);

% Define function handle to get the osculating state of chief and deputy
stateFuncD = @(tt) deval(stateHistDictD, tt);
stateFuncC = @(tt) deval(stateHistDictC, tt);


% Use numerical quadrature to compute mean state of the deputy and chief at each value of the argument of latitude
integralTol = 1e-10;
propagatedMeanStateHistC = zeros(length(ttHistory),6);
propagatedMeanStateHistD = zeros(length(ttHistory),6);
for i = 1:length(ttHistory)
    propagatedMeanStateHistD(i,:) = 1/(2*pi) * integral(stateFuncD, ttHistory(i)-pi, ttHistory(i)+pi, ...
        'RelTol', integralTol, 'AbsTol', integralTol, 'ArrayValued', true);
    propagatedMeanStateHistC(i,:) = 1/(2*pi) * integral(stateFuncC, ttHistory(i)-pi, ttHistory(i)+pi, ...
        'RelTol', integralTol, 'AbsTol', integralTol, 'ArrayValued', true);
end

% Get mean relative state history
propagatedMeanRelStateHist = propagatedMeanStateHistD - propagatedMeanStateHistC;

%% Compute osculating state by applying the transformation mean->osculating to numerical propagation

% We take the numerical mean relative and absolute state histories: relMeanState(tt) and meanStateC(tt)
% Then for each value of tt, we transform: [relMeanState(tt), meanStateC(tt), tt] -> relOscState(tt)

% Define analytical state histories
mean2oscAnalyticalRelStateHist = zeros(length(ttHistory),6); % State history based on full analytical solution

for i = 1:length(ttHistory)
    
    % Apply analytical transformation from mean to osculating elements
    mean2oscAnalyticalRelStateHist(i,:) = meanToOscRelativeState(...
        ttHistory(i), mu, j2, R, propagatedMeanRelStateHist(i,:), propagatedMeanStateHistC(i,:), kD, rhoD, kC, rhoC, expansionOrder);
end

%% Plot numerically propagated osculating state history, and analytical mean state history

% Define labels, units, and scaling factors for plotting
y_label = ["\delta \beta", "\delta e_x", "\delta e_y", "\delta p", "\delta \Omega", "\delta t"];
unit_label = [" [-]", " [-]", " [-]", " [-]", " [deg]", " [s]"];
factor = [1, j2, j2, 1, 180/pi, 1];
num_vars = length(y_label);

% Convert argument of latitude to revolutions
ttHistoryRev = ttHistory/(2*pi);

% Plot font size
font_size = 14;

try
    close(1)
catch
end
figure(1);
set(gcf,'Position', [100 100 800 425])
for var = 1:num_vars
    subplot(2,3,var)

    hold on
    plot(ttHistoryRev, propagatedOsculatingRelStateHist(:,var)*factor(var), 'linewidth', 2, "DisplayName", "Osc (numerical)");
    plot(ttHistoryRev, osc2meanAnalyticalRelStateHist(:,var)*factor(var), 'linewidth', 2, "DisplayName", "Osc2mean (analytical)");
    hold off

    grid;
    xlabel('\theta [2\pi]')
    ylabel(y_label(var) + unit_label(var))
    set(gca, 'FontSize', font_size);

    xlim([min(ttHistoryRev), max(ttHistoryRev)])
    
    if var == 4
        legend("Location", "best")
    end
end

sgtitle("Osculating to mean transformation")

%% Plot numerical mean state history, and analytical osculating state history

try
    close(2)
catch
end
figure(2);
set(gcf,'Position', [100 100 800 425])
for var = 1:num_vars
    subplot(2,3,var)

    hold on
    plot(ttHistoryRev, propagatedMeanRelStateHist(:,var)*factor(var), 'linewidth', 2, "DisplayName", "Mean (numerical)");
    plot(ttHistoryRev, mean2oscAnalyticalRelStateHist(:,var)*factor(var), 'linewidth', 2, "DisplayName", "Mean2osc (analytical)");
    hold off

    grid;
    xlabel('\theta [2\pi]')
    ylabel(y_label(var) + unit_label(var))
    set(gca, 'FontSize', font_size);

    xlim([min(ttHistoryRev), max(ttHistoryRev)])
    
    if var == 4
        legend("Location", "best")
    end
end

sgtitle("Mean to osculating transformation")