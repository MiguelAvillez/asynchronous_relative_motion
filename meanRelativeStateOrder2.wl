(* ::Package:: *)

(*
This script uses the nomenclature, formulations and solutions from: 
M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion 
Approximations Under J2 and Atmospheric Drag", TODO

Summary: Provides the 2nd order relative state mean solution.

Notation:
	tt: argument of latitude at which to compute the mean relative state
    mu: gravitational parameter
    j2: J2 coefficient of the gravity model
    R: Radius of the central planet
    tti: Initial argument of latitude
    kD: Ballistic parameter of the deputy (D)
    rhoD: Atmospheric density acting on the deputy (D)
    kC: Ballistic parameter of the chief (C)
    rhoC: Atmospheric density acting on the chief (C)
    [dbbMi; dxMi; dyMi; dpMi; dooMi; dtMi]: initial mean relative state
    [bbMi; xMi; yMi; pMi; ooMi; tMi]: initial mean state of the deputy

Solution can be imported to a Mathematica notebook (.nb extension) using 
the Import[...] command.
*)

dbbM[tt_] = dbbMi - (R*(tt - tti)*(4*bbMi*kC*rhoC - 4*bbMi*kD*rhoD + 
        4*dbbMi*kD*rhoD + 3*bbMi^5*j2*(-1 + bbMi^2*pMi^2)*
         (-(kC*rhoC) + kD*rhoD)*Cos[2*tti]))/(8*bbMi^2)
 
dxM[tt_] = dxMi + (3*bbMi^2*(tt - tti)*
       (bbMi*j2*(4*dbbMi*yMi + bbMi*(dyMi - 5*bbMi^2*dyMi*pMi^2 - 
            10*bbMi*pMi*(bbMi*dpMi + 3*dbbMi*pMi)*yMi)) - 
        (-1 + 5*bbMi^2*pMi^2)*R*(-(kC*rhoC) + kD*rhoD)*Cos[tti]))/4
 
dyM[tt_] = dyMi + (3*bbMi^2*(tt - tti)*
       (bbMi*j2*(-(bbMi*dxMi) - 4*dbbMi*xMi + 30*bbMi^2*dbbMi*pMi^2*xMi + 
          5*bbMi^3*pMi*(dxMi*pMi + 2*dpMi*xMi)) - (-1 + 5*bbMi^2*pMi^2)*R*
         (-(kC*rhoC) + kD*rhoD)*Sin[tti]))/4
 
dpM[tt_] = dpMi + (R*(tt - tti)*(2*bbMi*kC*pMi*rhoC + 4*dbbMi*kD*pMi*rhoD - 
        2*bbMi*kD*(dpMi + pMi)*rhoD + 3*bbMi^5*j2*pMi*(-1 + bbMi^2*pMi^2)*
         (-(kC*rhoC) + kD*rhoD)*Cos[2*tti]))/(4*bbMi^3)
 
dooM[tt_] = dooMi + (3*bbMi^3*j2*(-2*bbMi*dpMi*(2*bbMi + 10*dbbMi + 
          bbMi^5*j2) - 2*dbbMi*(10*bbMi + 20*dbbMi + 9*bbMi^5*j2)*pMi + 
        15*bbMi^8*dpMi*j2*pMi^2 + 55*bbMi^7*dbbMi*j2*pMi^3)*(tt - tti))/8
 
dtM[tt_] = (128*bbMi^7*dtMi*Sqrt[mu] + 135*bbMi^10*j2^2*kC*Pi^2*
       (-1 + bbMi^2*pMi^2)^2*R^(5/2)*rhoC + 160*bbMi*dbbMi*kD*Pi^2*R^(5/2)*
       rhoD - 480*dbbMi^2*kD*Pi^2*R^(5/2)*rhoD - 135*bbMi^10*j2^2*kD*Pi^2*
       R^(5/2)*rhoD + 270*bbMi^12*j2^2*kD*Pi^2*pMi^2*R^(5/2)*rhoD - 
      135*bbMi^14*j2^2*kD*Pi^2*pMi^4*R^(5/2)*rhoD + 120*bbMi^6*j2^2*Pi^2*
       (-1 + bbMi^2*pMi^2)*R^(5/2)*(-(kC*rhoC) + kD*rhoD)*xMi*Cos[tti] + 
      192*bbMi^6*j2*R^(3/2)*Cos[tt]*(bbMi*j2*(-(bbMi*dyMi) - dbbMi*yMi + 
          15*bbMi^2*dbbMi*pMi^2*yMi + 5*bbMi^3*pMi*(dyMi*pMi + 2*dpMi*yMi)) + 
        (-1 + 5*bbMi^2*pMi^2)*R*(-(kC*rhoC) + kD*rhoD)*Cos[tti]) + 
      2*bbMi^5*j2*R^(3/2)*(15*Pi^2*R*(-(bbMi*kC*(-1 + bbMi^2*pMi^2)*
            (4 - 3*bbMi^4*j2 + 3*bbMi^6*j2*pMi^2)*rhoC) + 
          kD*(-4*bbMi + 4*dbbMi + 3*bbMi^5*j2 + 8*bbMi^3*dpMi*pMi + 
            2*bbMi^2*(2*bbMi + 2*dbbMi - 3*bbMi^5*j2)*pMi^2 + 
            3*bbMi^9*j2*pMi^4)*rhoD)*Cos[2*tti] + 96*bbMi*Sin[tt]*
         (bbMi*j2*(dbbMi*xMi + bbMi*(dxMi - 5*bbMi^2*dxMi*pMi^2 - 
              5*bbMi*pMi*(2*bbMi*dpMi + 3*dbbMi*pMi)*xMi)) + 
          (-1 + 5*bbMi^2*pMi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[tti]) + 
        5*bbMi*j2*Pi^2*(-1 + bbMi^2*pMi^2)*R*(-(kC*rhoC) + kD*rhoD)*
         (4*xMi*Cos[3*tti] + 3*bbMi^4*(2 + bbMi^2*pMi^2)*Cos[4*tti] - 
          16*yMi*Sin[tti]^3)) - 3*R^(5/2)*(tt - tti)^2*
       (-(bbMi^2*kC*(32 + 135*bbMi^8*j2^2*(-1 + bbMi^2*pMi^2)^2)*rhoC) + 
        kD*(32*bbMi^2 - 160*bbMi*dbbMi + 480*dbbMi^2 + 135*bbMi^10*j2^2 - 
          270*bbMi^12*j2^2*pMi^2 + 135*bbMi^14*j2^2*pMi^4)*rhoD - 
        10*bbMi^5*j2*(12*bbMi*j2*(-1 + bbMi^2*pMi^2)*(-(kC*rhoC) + kD*rhoD)*
           xMi*Cos[tti] + 3*(-(bbMi*kC*(-1 + bbMi^2*pMi^2)*(4 - 3*bbMi^4*j2 + 
               3*bbMi^6*j2*pMi^2)*rhoC) + kD*(-4*bbMi + 4*dbbMi + 
              3*bbMi^5*j2 + 8*bbMi^3*dpMi*pMi + 2*bbMi^2*(2*bbMi + 2*dbbMi - 
                3*bbMi^5*j2)*pMi^2 + 3*bbMi^9*j2*pMi^4)*rhoD)*Cos[2*tti] + 
          bbMi*j2*(-1 + bbMi^2*pMi^2)*(-(kC*rhoC) + kD*rhoD)*
           (4*xMi*Cos[3*tti] + 3*bbMi^4*(2 + bbMi^2*pMi^2)*Cos[4*tti] - 
            16*yMi*Sin[tti]^3))) + 4*bbMi*R^(3/2)*(tt - tti)*
       (192*bbMi*dbbMi^2 - 320*dbbMi^3 + 48*bbMi^6*dbbMi*j2 + 
        405*bbMi^10*dbbMi*j2^2 - 732*bbMi^13*dpMi*j2^2*pMi - 
        576*bbMi^7*dbbMi^2*j2*pMi^2 - 2562*bbMi^12*dbbMi*j2^2*pMi^2 + 
        1572*bbMi^15*dpMi*j2^2*pMi^3 + 3537*bbMi^14*dbbMi*j2^2*pMi^4 - 
        576*bbMi^8*dbbMi*j2*pMi*(2*dpMi + pMi) - 192*bbMi^9*dpMi*j2*
         (dpMi + 2*pMi) + 96*bbMi^3*j2^2*(dxMi*xMi + dyMi*yMi) - 
        48*bbMi^2*dbbMi*(2 + 3*j2^2*(xMi^2 + yMi^2)) + 
        3*bbMi*j2*(8*(3*bbMi^5*j2*(-(bbMi*dxMi) - 5*dbbMi*xMi + 
              7*bbMi^2*dbbMi*pMi^2*xMi + bbMi^3*pMi*(dxMi*pMi + 
                2*dpMi*xMi)) + 4*R*(-(kC*rhoC) + kD*rhoD)*yMi)*Cos[tti] + 
          6*bbMi^4*(dbbMi^2*(-40 + 84*bbMi^2*pMi^2) + 
            dbbMi*(27*bbMi^5*j2 + 56*bbMi^3*dpMi*pMi - 66*bbMi^7*j2*pMi^2 + 
              39*bbMi^9*j2*pMi^4) + 4*bbMi^4*dpMi*(dpMi + 3*bbMi^4*j2*pMi*(
                -1 + bbMi^2*pMi^2)))*Cos[2*tti] + 8*bbMi^5*j2*
           (-(bbMi*dxMi) - 5*dbbMi*xMi + 7*bbMi^2*dbbMi*pMi^2*xMi + 
            bbMi^3*pMi*(dxMi*pMi + 2*dpMi*xMi))*Cos[3*tti] + 
          3*bbMi^9*j2*(4*bbMi^3*dpMi*pMi*(-8 + 11*bbMi^2*pMi^2) + 
            dbbMi*(45 - 176*bbMi^2*pMi^2 + 143*bbMi^4*pMi^4))*Cos[4*tti] - 
          8*(4*R*(-(kC*rhoC) + kD*rhoD)*xMi + 3*bbMi^5*j2*(-(bbMi*dyMi) - 
              5*dbbMi*yMi + 7*bbMi^2*dbbMi*pMi^2*yMi + bbMi^3*pMi*(dyMi*pMi + 
                2*dpMi*yMi)))*Sin[tti] - 16*bbMi^4*(-1 + bbMi^2*pMi^2)*R*
           (-(kC*rhoC) + kD*rhoD)*Sin[2*tti] + 8*bbMi^5*j2*
           (-(bbMi*dyMi) - 5*dbbMi*yMi + 7*bbMi^2*dbbMi*pMi^2*yMi + 
            bbMi^3*pMi*(dyMi*pMi + 2*dpMi*yMi))*Sin[3*tti])))/
     (128*bbMi^7*Sqrt[mu])
