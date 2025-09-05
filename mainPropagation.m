%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion Approximations 
%   Under J2 and Atmospheric Drag", TODO
%
% Summary: Applies the solutions from section IV, V.C, and VI of the paper
%   Propagates the relative state using:
%       - Full analytical osculating-elements solution
%       - STM-based osculating-elements solution
%       - Full analytical mean-elements solution
%       - STM-based mean-elements solution
%   Plots the computed state evolutions
%
% Dependencies:
%   osculatingRelativeState.m
%   osculatingRelativeStateStm.m
%   osculatingRelativeAltitude.m
%   meanRelativeState.m
%   meanRelativeStateStm.m
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
%   Plot of relative orbital elements (osculating)
%   Plot of relative orbital elements error (osculating)
%   Plot of relative altitude vs along-track and cross-track (osculating)
%   Plot of relative orbital elements (mean)
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
numRevolutions = 200;
if numRevolutions <= 10
    numPointsPerRev = 100;
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

%% Numerically propagate the state of the chief and deputy

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

% Get relative state history
propagatedRelStateHist = propagatedStateHistD - propagatedStateHistC;

%% Compute relative state and relative altitude using analytical osculating solutions

% Define analytical state histories
osculatingAnalyticalRelStateHist = zeros(length(ttHistory),6); % State history based on full analytical solution
osculatingStmRelStateHist = zeros(length(ttHistory),6); % State history based on osculating STM

% Define analytical relative altitude history
osculatingAnalyticalDrHist = zeros(size(ttHistory));

for i = 1:length(ttHistory)
    % Compute osculating state based on full analytical solution
    osculatingAnalyticalRelStateHist(i,:) = osculatingRelativeState(...
        ttHistory(i), mu, j2, R, tti, relativeStateInitial, stateInitialChief, kD, rhoD, kC, rhoC, expansionOrder);

    % Compute osculating state based on STM solution
    % First get STM and sensitivity matrix, then multiply them by the deviations
    [phi, s] = osculatingRelativeStateStm(ttHistory(i), mu, j2, R, tti, stateInitialChief, kC, rhoC, expansionOrder);
    osculatingStmRelStateHist(i,:) = phi * relativeStateInitial + s * dzeta;

    % Compute the relative altitude
    osculatingAnalyticalDrHist(i) = osculatingRelativeAltitude(...
        ttHistory(i), mu, j2, R, tti, relativeStateInitial, stateInitialChief, kD, rhoD, kC, rhoC, expansionOrder);

end

%% Numerically compute initial state in mean elements

% Numerically propagate full dynamics: backward for pi
[~, auxStateHistD] = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionD),...
    [ttHistory(1), ttHistory(1)-pi], stateInitialDeputy, odeOptions);
[~, auxStateHistC] = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionC),...
    [ttHistory(1), ttHistory(1)-pi], stateInitialChief, odeOptions);

% Propagate dynamics over range [-pi, pi], which is range necessary to evaluate initial mean state
stateHistDictD = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionD),...
    [ttHistory(1)-pi, ttHistory(1)+pi], auxStateHistD(end,:), odeOptions);
stateHistDictC = ode89(@(t,x) stateArgLatDerivative(t, x, mu, R, j2, accFunctionC),...
    [ttHistory(1)-pi, ttHistory(1)+pi], auxStateHistC(end,:), odeOptions);

% Define function handle to integrate mean state
stateFuncD = @(tt) deval(stateHistDictD, tt);
stateFuncC = @(tt) deval(stateHistDictC, tt);

integralTol = 1e-10;
% Compute the initial mean state of the chief
stateInitialMeanChief = 1/(2*pi) * integral(stateFuncC, ttHistory(1)-pi, ttHistory(1)+pi, ...
    'RelTol', integralTol, 'AbsTol', integralTol, 'ArrayValued', true);
% Compute the initial mean state of the deputy
stateInitialMeanDeputy = 1/(2*pi) * integral(stateFuncD, ttHistory(1)-pi, ttHistory(1)+pi, ...
    'RelTol', integralTol, 'AbsTol', integralTol, 'ArrayValued', true);
% Compute the initial mean relative state
stateInitialMeanRelative = stateInitialMeanDeputy - stateInitialMeanChief;

%% Compute analytical mean elements solution

% Define analytical state histories
meanAnalyticalRelStateHist = zeros(length(ttHistory),6); % State history based on mean-elements full analytical solution
meanStmRelStateHist = zeros(length(ttHistory),6); % State history based on mean-elements STM solution

for i = 1:length(ttHistory)

    % Compute mean state based on full analytical solution
    meanAnalyticalRelStateHist(i,:) = meanRelativeState(...
        ttHistory(i), mu, j2, R, tti, stateInitialMeanRelative, stateInitialMeanChief, kD, rhoD, kC, rhoC, expansionOrder);
        
    % Compute mean state based on STM solution
    % First get STM and sensitivity matrix, then multiply them by the deviations
    [phi, s] = meanRelativeStateStm(ttHistory(i), mu, j2, R, tti, stateInitialMeanChief, kC, rhoC, expansionOrder);
    meanStmRelStateHist(i,:) = phi * stateInitialMeanRelative + s * dzeta;
end

%% Plot numerically propagated and analytical osculating state history

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
    plot(ttHistoryRev, propagatedRelStateHist(:,var)*factor(var), 'linewidth', 3, "DisplayName", "Numerical");
    plot(ttHistoryRev, osculatingAnalyticalRelStateHist(:,var)*factor(var), 'linewidth', 2, "LineStyle", "--", ...
        "DisplayName", "Analytical");

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

sgtitle("Osculating elements")

%% Plot state error from full analytical and STM osculating solutions

try
    close(2)
catch
end
figure(2)
set(gcf,'Position', [100 100 800 425])
for var = 1:num_vars

    subplot(2,3,var)

    hold on
    pC = plot(ttHistoryRev, (osculatingAnalyticalRelStateHist(:,var) - propagatedRelStateHist(:,var))*factor(var), ...
        "LineWidth", 2, "LineStyle", "--", "DisplayName", "Analytical Full");
    plot(ttHistoryRev, (osculatingStmRelStateHist(:,var) - propagatedRelStateHist(:,var))*factor(var), ...
        "LineWidth", 2, "DisplayName", "STM");
    hold off

    uistack(pC,'top');

    grid;
    xlabel('\theta [2\pi]')
    ylabel("Error " + y_label(var) + unit_label(var))
    set(gca, 'FontSize', font_size);

    xlim([min(ttHistoryRev), max(ttHistoryRev)])

    if var == 4
        legend("Location", "best")
    end
end

sgtitle("Osculating elements error");

%% Plot along-track, cross-track, relative altitude based on full analytical solution

try
    close(3)
catch
end
figure(3)
set(gcf,'Position', [100 100 800 325])

% Along-track vs relative altitude
subplot(1,2,1)
plot(osculatingAnalyticalRelStateHist(:,6), osculatingAnalyticalDrHist, "LineWidth", 2);
xlabel('\delta t [s]');

% Cross-track vs relative altitude
subplot(1,2,2)
plot((osculatingAnalyticalRelStateHist(:,5) - we*osculatingAnalyticalRelStateHist(:,6))*rad2deg, osculatingAnalyticalDrHist, "LineWidth", 2);
xlabel('\delta \lambda [deg]');

for i = 1:2
    subplot(1,2,i)
    grid;
    ylabel('\delta r [m]');
    set(gca, 'FontSize', font_size);
end

sgtitle("Along-track, cross-track, relative altitude")

%% Plot mean state history based on full analytical and STM solutions

RGB = get(groot, "FactoryAxesColorOrder");

try
    close(4)
catch
end
figure(4);
set(gcf,'Position', [100 100 800 425])

for var = 1:num_vars
    subplot(2,3,var)
    hold on

    plot(ttHistoryRev, osculatingAnalyticalRelStateHist(:,var)*factor(var), 'linewidth', 2, "Color", [0.7, 0.7, 0.7], ...
        "DisplayName", "Osculating (analytical)");
    plot(ttHistoryRev, meanAnalyticalRelStateHist(:,var)*factor(var), 'linewidth', 2, "LineStyle", "-", "Color", RGB(1,:), ...
        "DisplayName", "Mean (analytical)");
    plot(ttHistoryRev, meanStmRelStateHist(:,var)*factor(var), 'linewidth', 2, "LineStyle", "--", "Color", RGB(2,:), ...
        "DisplayName", "Mean (STM)");

    hold off
    grid;
    xlabel('\theta [2\pi]')
    ylabel(y_label(var) + unit_label(var))
    set(gca, 'FontSize', font_size);

    xlim([min(ttHistoryRev), max(ttHistoryRev)])

    if var == 4
        legend("Location", "south")
    end
end
sgtitle("Mean elements")