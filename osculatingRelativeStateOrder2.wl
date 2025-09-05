(* ::Package:: *)

(*
This script uses the nomenclature, formulations and solutions from: 
M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion 
Approximations Under J2 and Atmospheric Drag", TODO

Summary: Provides the 2nd order relative state osculating solution.

Notation:
	tt: argument of latitude at which to compute the osculating relative state
    mu: gravitational parameter
    j2: J2 coefficient of the gravity model
    R: Radius of the central planet
    tti: Initial argument of latitude
    kD: Ballistic parameter of the deputy (D)
    rhoD: Atmospheric density acting on the deputy (D)
    kC: Ballistic parameter of the chief (C)
    rhoC: Atmospheric density acting on the chief (C)
    [dbbi; dxi; dyi; dpi; dooi; dti]: initial osculating relative state
    [bbCi; xCi; yCi; pCi; ooCi; tCi]: initial osculating state of the deputy

Solution can be imported to a Mathematica notebook (.nb extension) using 
the Import[...] command.
*)

dbb[tt_] = dbbi + (R*(-(bbCi*kC*rhoC) + (bbCi - dbbi)*kD*rhoD)*(tt - tti))/
      (2*bbCi^2) + (3*bbCi^4*j2*(2*bbCi^3*dpi*pCi + 
        dbbi*(-5 + 7*bbCi^2*pCi^2))*(Cos[2*tt] - Cos[2*tti]))/4
 
dx[tt_] = dxi + (bbCi^3*(3*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 15*bbCi^2*pCi^2))*
         Cos[tt] - 7*(bbCi^3*dpi*pCi + dbbi*(-2 + 3*bbCi^2*pCi^2))*
         Cos[3*tt] + 2*Cos[tti]*(10*dbbi - 11*bbCi^3*dpi*pCi - 
          33*bbCi^2*dbbi*pCi^2 + 7*(-2*dbbi + bbCi^3*dpi*pCi + 
            3*bbCi^2*dbbi*pCi^2)*Cos[2*tti])))/4 - 
     (R*(-(kC*rhoC) + kD*rhoD)*(Sin[tt] - Sin[tti]))/(bbCi^2*j2) + 
     (3*bbCi^3*j2*(tt - tti)*(8*dbbi*yCi + 2*bbCi*(dyi - 5*bbCi^2*dyi*pCi^2 - 
          10*bbCi*pCi*(bbCi*dpi + 3*dbbi*pCi)*yCi) - 
        bbCi^4*(-3*(bbCi^3*dpi*pCi*(-11 + 35*bbCi^2*pCi^2) + 
            dbbi*(6 - 55*bbCi^2*pCi^2 + 105*bbCi^4*pCi^4))*Sin[tti] + 
          7*(bbCi^3*dpi*pCi*(-3 + 5*bbCi^2*pCi^2) + 
            dbbi*(2 - 15*bbCi^2*pCi^2 + 15*bbCi^4*pCi^4))*Sin[3*tti])))/8
 
dy[tt_] = dyi + (R*(-(kC*rhoC) + kD*rhoD)*(Cos[tt] - Cos[tti]))/(bbCi^2*j2) + 
     (3*bbCi^3*j2*(tt - tti)*(-2*bbCi*dxi + 10*bbCi^3*dxi*pCi^2 - 
        8*dbbi*xCi + 20*bbCi^3*dpi*pCi*xCi + 60*bbCi^2*dbbi*pCi^2*xCi - 
        3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + 
          dbbi*(-2 + 15*bbCi^2*pCi^2))*Cos[tti] + 
        7*bbCi^4*(bbCi^3*dpi*pCi*(-3 + 5*bbCi^2*pCi^2) + 
          dbbi*(2 - 15*bbCi^2*pCi^2 + 15*bbCi^4*pCi^4))*Cos[3*tti]))/8 + 
     bbCi^3*(6*dbbi*Sin[tt] + 7*(bbCi^3*dpi*pCi + dbbi*(-2 + 3*bbCi^2*pCi^2))*
        Sin[tt]^3 - 6*dbbi*Sin[tti] - 7*(-2*dbbi + bbCi^3*dpi*pCi + 
         3*bbCi^2*dbbi*pCi^2)*Sin[tti]^3)
 
dp[tt_] = dpi + ((2*dbbi*kD*pCi*R*rhoD + bbCi*R*(kC*pCi*rhoC - 
          kD*(dpi + pCi)*rhoD))*(tt - tti))/(2*bbCi^3)
 
doo[tt_] = dooi + (3*bbCi^3*j2*(tt - tti)*
       (-2*bbCi*dpi*(2*bbCi + 10*dbbi + bbCi^5*j2) - 
        2*dbbi*(10*bbCi + 20*dbbi + 9*bbCi^5*j2)*pCi + 
        15*bbCi^8*dpi*j2*pCi^2 + 55*bbCi^7*dbbi*j2*pCi^3 + 
        15*bbCi^5*j2*(-(bbCi*dpi) - 9*dbbi*pCi + 3*bbCi^3*dpi*pCi^2 + 
          11*bbCi^2*dbbi*pCi^3)*Cos[2*tti]))/8 + 
     (3*bbCi^4*j2*(bbCi*dpi + 5*dbbi*pCi)*(Sin[2*tt] - Sin[2*tti]))/4
 
dt[tt_] = (-48*R^(5/2)*(-(bbCi^2*kC*rhoC) + (bbCi^2 - 5*bbCi*dbbi + 
          15*dbbi^2)*kD*rhoD)*(tt - tti)^2 + 
      16*bbCi^2*(4*bbCi^5*dti*Sqrt[mu] + 8*R^(5/2)*(-(kC*rhoC) + kD*rhoD) + 
        R^(3/2)*(bbCi*j2*(-8*(bbCi*dyi - 3*dbbi*yCi)*Cos[tti] - 
            (8*bbCi*dxi - 24*dbbi*xCi + bbCi^4*(3*dbbi*Cos[tti] - 
                bbCi^2*pCi*(14*bbCi*dpi*Cos[tt] + 15*(2*bbCi*dpi + 
                    3*dbbi*pCi)*Cos[tti]) + 7*(2*bbCi^3*dpi*pCi + 
                  dbbi*(-1 + 3*bbCi^2*pCi^2))*Cos[3*tti]))*Sin[tt]) + 
          (bbCi*j2*(8*(bbCi*dxi - 3*dbbi*xCi) + (9*bbCi^4*dbbi - 
                30*bbCi^7*dpi*pCi - 45*bbCi^6*dbbi*pCi^2)*Cos[tti]) + 
            8*R*(kC*rhoC - kD*rhoD)*Sin[tt])*Sin[tti]) + 
        R^(3/2)*Cos[tt]*(8*bbCi*j2*(bbCi*dyi - 3*dbbi*yCi) + 
          8*R*(kC*rhoC - kD*rhoD)*Cos[tti] - bbCi^5*j2*
           ((dbbi - 21*bbCi^2*dbbi*pCi^2)*Sin[tt] - 
            2*(dbbi - 14*bbCi^3*dpi*pCi - 21*bbCi^2*dbbi*pCi^2 + 
              7*(-dbbi + 2*bbCi^3*dpi*pCi + 3*bbCi^2*dbbi*pCi^2)*Cos[2*tti])*
             Sin[tti]))) + bbCi*R^(3/2)*(tt - tti)*(384*bbCi*dbbi^2 - 
        640*dbbi^3 + 96*bbCi^6*dbbi*j2 + 1185*bbCi^10*dbbi*j2^2 - 
        1980*bbCi^13*dpi*j2^2*pCi - 1152*bbCi^7*dbbi^2*j2*pCi^2 - 
        6930*bbCi^12*dbbi*j2^2*pCi^2 + 4740*bbCi^15*dpi*j2^2*pCi^3 + 
        10665*bbCi^14*dbbi*j2^2*pCi^4 - 1152*bbCi^8*dbbi*j2*pCi*
         (2*dpi + pCi) - 384*bbCi^9*dpi*j2*(dpi + 2*pCi) + 
        192*bbCi^3*j2^2*(dxi*xCi + dyi*yCi) - 96*bbCi^2*dbbi*
         (2 + 3*j2^2*(xCi^2 + yCi^2)) - 3*bbCi^6*j2*
         (24*j2*(bbCi*dxi + dbbi*xCi + 9*bbCi^2*dbbi*pCi^2*xCi + 
            3*bbCi^3*pCi*(dxi*pCi + 2*dpi*xCi))*Cos[tti] + 
          48*(dbbi - bbCi^3*dpi^2 + 5*bbCi^4*dbbi*j2 - 2*bbCi^2*dpi*
             (bbCi + 3*dbbi + bbCi^5*j2)*pCi - bbCi*dbbi*(3*bbCi + 3*dbbi + 
              7*bbCi^5*j2)*pCi^2)*Cos[2*tti] + 72*j2*(bbCi*dxi + dbbi*xCi - 
            3*bbCi^2*dbbi*pCi^2*xCi - bbCi^3*pCi*(dxi*pCi + 2*dpi*xCi))*
           Cos[3*tti] + 4*j2*Cos[tt]*(8*dbbi*xCi + 
            8*bbCi*(dxi - 5*bbCi^2*dxi*pCi^2 - 5*bbCi*pCi*(2*bbCi*dpi + 
                3*dbbi*pCi)*xCi) + 15*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
             (4*bbCi^3*dpi*pCi + dbbi*(-1 + 9*bbCi^2*pCi^2))*Cos[tti] + 
            7*bbCi^4*(4*bbCi^3*dpi*pCi*(3 - 5*bbCi^2*pCi^2) + 
              dbbi*(-5 + 42*bbCi^2*pCi^2 - 45*bbCi^4*pCi^4))*Cos[3*tti]) - 
          45*bbCi^4*j2*(-1 + bbCi^2*pCi^2)*(4*bbCi^3*dpi*pCi + 
            dbbi*(-5 + 9*bbCi^2*pCi^2))*Cos[4*tti] + 
          32*j2*(bbCi*dyi + dbbi*yCi - 15*bbCi^2*dbbi*pCi^2*yCi - 
            5*bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi))*Sin[tt] - 
          12*j2*(2*(5*bbCi*dyi + 5*dbbi*yCi - 27*bbCi^2*dbbi*pCi^2*yCi - 
              9*bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi)) + 
            bbCi^4*(4*bbCi^3*dpi*pCi*(11 - 35*bbCi^2*pCi^2) + 
              dbbi*(-15 + 154*bbCi^2*pCi^2 - 315*bbCi^4*pCi^4))*Sin[tt])*
           Sin[tti] + 4*j2*(18*(bbCi*dyi + dbbi*yCi - 3*bbCi^2*dbbi*pCi^
                2*yCi - bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi)) + 
            7*bbCi^4*(4*bbCi^3*dpi*pCi*(3 - 5*bbCi^2*pCi^2) + 
              dbbi*(-5 + 42*bbCi^2*pCi^2 - 45*bbCi^4*pCi^4))*Sin[tt])*
           Sin[3*tti])))/(64*bbCi^7*Sqrt[mu])
