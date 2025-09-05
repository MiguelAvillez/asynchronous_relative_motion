This repository contains the source-code describing the formulations and solutions in
	M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion Approximations 
	Under J2 and Atmospheric Drag", TODO

The files listed below give all analytical solutions proposed in the paper (for order 2 and 3):
1) Evolution of osculating relative orbital elements (as MATLAB and Mathematica code) 
2) Evolution of mean relative orbital elements (as MATLAB and Mathematica code) 
3) Evolution of osculating relative altitude (as MATLAB and Mathematica code) 
4) Evolution of osculating relative orbital elements based on STM (as MATLAB code) 
5) Evolution of mean relative orbital elements based on STM (as MATLAB code) 
6) Transformation from osculating to mean relative orbital elements (as MATLAB code) 
7) Transformation from mean to osculating relative orbital elements (as MATLAB code) 

The MATLAB code is applied in mainPropagation.m and mainTransformations.m

-----------------------------------------------------------------------------------
Nomenclature:
	Independent variable: tt
	Osculating absolute states: [bb; x; y; p; oo; t]
	Osculating relative states: [dbb; dx; dy; dp; doo; dt]
	Mean absolute states: [bbM; xM; yM; pM; ooM; tM]
	Mean relative states: [dbbM; dxM; dyM; dpM; dooM; dtM]

	Correspondence to notation used in the paper (LaTeX):
		tt -> \theta

		bb -> \beta
		x -> X
		y -> Y
		p -> p
		oo -> \Omega
		t -> t

		dbb -> \delta \beta
		dx -> \delta X
		dy -> \delta Y
		dp -> \delta p
		doo -> \delta \Omega
		dt -> \delta t

		bbM -> \bar{\beta}
		xM -> \bar{X}
		yM -> \bar{Y}
		pM -> \bar{p}
		ooM -> \bar{\Omega}
		tM -> \bar{t}

		dbbM -> \bar{\delta \beta}
		dxM -> \bar{\delta X}
		dyM -> \bar{\delta Y}
		dpM -> \bar{\delta p}
		dooM -> \bar{\delta \Omega}
		dtM -> \bar{\delta t}

-----------------------------------------------------------------------------------
Main files:
- mainPropagation.m
	- Applies the solutions from sections IV, V.C, and VI of the paper, that is:
		1) Propagation of relative osculating elements based on full analytical solution
		2) Propagation of relative osculating altitude based on full analytical solution
		3) Propagation of relative osculating elements based on STM solution
		4) Propagation of relative mean elements based on full analytical solution
		5) Propagation of relative mean elements based on STM solution

- mainTransformations.m
	- Applies the solutions from section V.A and V.B of the paper, that is:
		1) Analytical transformation from osculating to mean relative orbital elements
		2) Analytical transformation from mean to osculating relative orbital elements

-----------------------------------------------------------------------------------
Full list of files:

- mainPropagation.m: see above
- mainTransformations.m: see above

- accAtmosphericDrag.m: computes acceleration due to atmospheric drag
- accGravitationalJ2.m: computes acceleration due to J2
- stateArgLatDerivative.m: computes state derivative

- meanToOscAbsoluteState.m: transforms mean to osculating absolute states, for order 2 and 3
- meanToOscRelativeState.m: transforms mean to osculating relative states, for order 2 and 3
- oscToMeanAbsoluteState.m: transforms osculating to mean absolute states, for order 2 and 3
- oscToMeanRelativeState.m: transforms osculating to mean relative states, for order 2 and 3

- meanRelativeState.m: propagates mean relative state (full analytical solution), for order 2 and 3
- meanRelativeStateStm.m: propagates mean relative state (STM solution), for order 2 and 3
- osculatingRelativeAltitude.m: propagates osculating relative altitude (full analytical solution), for order 2 and 3
- osculatingRelativeState.m: propagates osculating relative state (full analytical solution), for order 2 and 3
- osculatingRelativeStateStm.m: propagates osculating relative state (STM solution), for order 2 and 3

- osculatingRelativeStateOrder2.wl: 2nd order osculating relative state solution (Mathematica code)
- osculatingRelativeStateOrder3.wl: 3rd order osculating relative state solution (Mathematica code)
- osculatingRelativeAltitudeOrder2.wl: 2nd order osculating relative altitude solution (Mathematica code)
- osculatingRelativeAltitudeOrder3.wl: 3rd order osculating relative altitude solution (Mathematica code)
- meanRelativeStateOrder2.wl: 2nd order mean relative state solution (Mathematica code)
- meanRelativeStateOrder3.wl: 3rd order mean relative state solution (Mathematica code)

