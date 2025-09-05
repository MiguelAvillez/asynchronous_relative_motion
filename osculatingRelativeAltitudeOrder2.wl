(* ::Package:: *)

(*
This script uses the nomenclature, formulations and solutions from: 
M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion 
Approximations Under J2 and Atmospheric Drag", TODO

Summary: Provides the 2nd order relative altitude osculating solution.

Notation:
	tt: argument of latitude at which to compute the osculating relative altitude
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

dr[tt_] = (R*(3*R*(-(kC*rhoC) + kD*rhoD)*(tt - tti)^2*
        (8*R*(kC*rhoC + kD*rhoD) + bbCi^6*j2^2*(-1 + 5*bbCi^2*pCi^2)*
          (-8*yCi*Cos[tt] + 8*xCi*Sin[tt] + bbCi^4*(3*(1 - 5*bbCi^2*pCi^2)*
              Cos[tti]*Sin[tt] + 7*(-1 + bbCi^2*pCi^2)*Sin[tt - 3*tti] + 
             3*(-3 + 7*bbCi^2*pCi^2)*Cos[tt]*Sin[tti]))) + 
       2*bbCi*(tt - tti)*(-36*bbCi^5*j2*(-1 + bbCi^2*pCi^2)*R*
          (-(kC*rhoC) + kD*rhoD)*Cos[2*tti] - 3*bbCi^5*j2*
          (1 - 5*bbCi^2*pCi^2)*Cos[tti]*(2*R*(kC*rhoC - kD*rhoD)*Cos[tt] + 
           3*bbCi^5*j2*(10*bbCi^3*dpi*pCi + dbbi*(-3 + 25*bbCi^2*pCi^2))*
            Sin[tt]) + 7*bbCi^5*j2*Cos[3*tti]*(2*(-1 + bbCi^2*pCi^2)*R*
            (-(kC*rhoC) + kD*rhoD)*Cos[tt] - 3*bbCi^5*j2*
            (2*bbCi^3*dpi*pCi*(-3 + 5*bbCi^2*pCi^2) + 
             dbbi*(3 - 24*bbCi^2*pCi^2 + 25*bbCi^4*pCi^4))*Sin[tt]) + 
         4*(bbCi*kC*(4 + 3*bbCi^4*j2 - 9*bbCi^6*j2*pCi^2)*R*rhoC + 
           kD*(-4*bbCi + 16*dbbi - 3*bbCi^5*j2 + 9*bbCi^7*j2*pCi^2)*R*rhoD + 
           bbCi*j2*((4*R*(-(kC*rhoC) + kD*rhoD)*xCi + 3*bbCi^5*j2*
                (-(bbCi*dyi) - 2*dbbi*yCi + 20*bbCi^2*dbbi*pCi^2*yCi + 
                 5*bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi)))*Cos[tt] + 
             4*bbCi^4*(-1 + bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*
              Cos[2*tt] + (3*bbCi^5*j2*(2*dbbi*xCi + bbCi*(dxi - 5*bbCi^2*dxi*
                    pCi^2 - 10*bbCi*pCi*(bbCi*dpi + 2*dbbi*pCi)*xCi)) + 4*R*
                (-(kC*rhoC) + kD*rhoD)*yCi)*Sin[tt])) - 
         3*bbCi^5*j2*(3*bbCi^5*j2*(-22*bbCi^3*dpi*pCi + 70*bbCi^5*dpi*pCi^3 + 
             dbbi*(9 - 88*bbCi^2*pCi^2 + 175*bbCi^4*pCi^4))*Cos[tt] + 
           2*(-3 + 7*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[tt])*
          Sin[tti] + 7*bbCi^5*j2*(3*bbCi^5*j2*(2*bbCi^3*dpi*pCi*
              (-3 + 5*bbCi^2*pCi^2) + dbbi*(3 - 24*bbCi^2*pCi^2 + 25*bbCi^4*
                pCi^4))*Cos[tt] + 2*(-1 + bbCi^2*pCi^2)*R*(-(kC*rhoC) + 
             kD*rhoD)*Sin[tt])*Sin[3*tti]) + 
       8*bbCi^2*(2*bbCi^5*j2*(dbbi - bbCi^3*dpi*pCi - 2*bbCi^2*dbbi*pCi^2)*
          Cos[2*tt] + 2*(dbbi*(-4*bbCi + 6*dbbi + 3*bbCi^5*j2) - 
           9*bbCi^8*dpi*j2*pCi - 18*bbCi^7*dbbi*j2*pCi^2 + 
           6*bbCi^5*j2*(bbCi^3*dpi*pCi + dbbi*(-1 + 2*bbCi^2*pCi^2))*
            Cos[2*tti]) + Cos[tt]*(-4*bbCi*j2*(bbCi*dxi - 2*dbbi*xCi) + 
           bbCi^5*j2*(3*(5*bbCi^3*dpi*pCi + dbbi*(-1 + 10*bbCi^2*pCi^2))*
              Cos[tti] + 7*(dbbi - bbCi^3*dpi*pCi - 2*bbCi^2*dbbi*pCi^2)*
              Cos[3*tti]) + 4*R*(kC*rhoC - kD*rhoD)*Sin[tti]) + 
         Sin[tt]*(4*R*(-(kC*rhoC) + kD*rhoD)*Cos[tti] + 
           bbCi*j2*(-4*bbCi*dyi + 8*dbbi*yCi + 3*bbCi^4*(7*bbCi^3*dpi*pCi + 
               dbbi*(-3 + 14*bbCi^2*pCi^2))*Sin[tti] + 7*bbCi^4*
              (dbbi - bbCi^3*dpi*pCi - 2*bbCi^2*dbbi*pCi^2)*Sin[3*tti])))))/
     (32*bbCi^6)
