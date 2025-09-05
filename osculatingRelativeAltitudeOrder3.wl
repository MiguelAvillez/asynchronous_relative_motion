(* ::Package:: *)

(*
This script uses the nomenclature, formulations and solutions from: 
M. Avillez and D. Arnas, "Osculating and Mean Asynchronous Relative Motion 
Approximations Under J2 and Atmospheric Drag", TODO

Summary: Provides the 3rd order relative altitude osculating solution.

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

dr[tt_] = (R*(-243*bbCi^19*j2^6*(1 - 5*bbCi^2*pCi^2)^4*R*
        (-(kC*rhoC) + kD*rhoD)*(tt - tti)^5*
        (Cos[tt]*(8*xCi + 3*bbCi^4*(1 - 5*bbCi^2*pCi^2)*Cos[tti]) + 
          8*yCi*Sin[tt] + bbCi^4*(7*(-1 + bbCi^2*pCi^2)*Cos[tt - 3*tti] + 
            3*(3 - 7*bbCi^2*pCi^2)*Sin[tt]*Sin[tti]))^2 + 
       162*bbCi^9*j2^3*(1 - 5*bbCi^2*pCi^2)^2*(tt - tti)^4*
        (Cos[tt]*(8*xCi + 3*bbCi^4*(1 - 5*bbCi^2*pCi^2)*Cos[tti]) + 
         8*yCi*Sin[tt] + bbCi^4*(7*(-1 + bbCi^2*pCi^2)*Cos[tt - 3*tti] + 
           3*(3 - 7*bbCi^2*pCi^2)*Sin[tt]*Sin[tti]))*
        (32*R^2*(-(kC*rhoC) + kD*rhoD)*(kC*rhoC + kD*rhoD) + 
         bbCi^6*j2^2*(-1 + 5*bbCi^2*pCi^2)*
          (8*(3*bbCi^5*j2*(-(bbCi*dxi) - 7*dbbi*xCi + 55*bbCi^2*dbbi*pCi^2*
                xCi + 5*bbCi^3*pCi*(dxi*pCi + 4*dpi*xCi)) + 
             8*R*(kC*rhoC - kD*rhoD)*yCi)*Cos[tt] + 
           8*(8*R*(-(kC*rhoC) + kD*rhoD)*xCi + 3*bbCi^5*j2*(-(bbCi*dyi) - 7*
                dbbi*yCi + 55*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*
                (dyi*pCi + 4*dpi*yCi)))*Sin[tt] - 3*bbCi^4*
            (-1 + 5*bbCi^2*pCi^2)*Cos[tti]*(3*bbCi^5*j2*(30*bbCi^3*dpi*pCi + 
               dbbi*(-11 + 85*bbCi^2*pCi^2))*Cos[tt] + 
             8*R*(-(kC*rhoC) + kD*rhoD)*Sin[tt]) + 7*bbCi^4*Cos[3*tti]*
            (3*bbCi^5*j2*(-22*bbCi^3*dpi*pCi + 30*bbCi^5*dpi*pCi^3 + dbbi*
                (11 - 88*bbCi^2*pCi^2 + 85*bbCi^4*pCi^4))*Cos[tt] + 
             8*(-1 + bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[tt]) - 
           3*bbCi^4*(-8*(-3 + 7*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*
              Cos[tt] + 3*bbCi^5*j2*(-74*bbCi^3*dpi*pCi + 210*bbCi^5*dpi*
                pCi^3 + dbbi*(33 - 316*bbCi^2*pCi^2 + 595*bbCi^4*pCi^4))*
              Sin[tt])*Sin[tti] + 7*bbCi^4*(-8*(-1 + bbCi^2*pCi^2)*R*
              (-(kC*rhoC) + kD*rhoD)*Cos[tt] + 3*bbCi^5*j2*(-22*bbCi^3*dpi*
                pCi + 30*bbCi^5*dpi*pCi^3 + dbbi*(11 - 88*bbCi^2*pCi^2 + 
                 85*bbCi^4*pCi^4))*Sin[tt])*Sin[3*tti])) - 
       (6*(tt - tti)^3*(49152*kC^2*R^3*rhoC^2*(-(kC*rhoC) + kD*rhoD) + 
          16384*R^3*(-(kC*rhoC) + kD*rhoD)^3 - 2304*bbCi^11*j2^3*
           (1 - 5*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Cos[tt]*
           (4*(-(bbCi*dxi) - 8*dbbi*xCi + 60*bbCi^2*dbbi*pCi^2*xCi + 
              5*bbCi^3*pCi*(dxi*pCi + 4*dpi*xCi)) - 9*bbCi^4*
             (-1 + 5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 
                15*bbCi^2*pCi^2))*Cos[tti] + 7*bbCi^4*
             (bbCi^3*dpi*pCi*(-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 
                47*bbCi^2*pCi^2 + 45*bbCi^4*pCi^4))*Cos[3*tti]) + 
          1152*bbCi^12*j2*(j2 - 5*bbCi^2*j2*pCi^2)^2*R*(-(kC*rhoC) + kD*rhoD)*
           (7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[tt - 3*tti] + 
            Cos[tt]*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti]) + 
            Sin[tt]*(8*yCi + 9*bbCi^4*Sin[tti] - 21*bbCi^6*pCi^2*Sin[tti])) - 
          405*bbCi^20*j2^6*(1 - 5*bbCi^2*pCi^2)^4*R*(-(kC*rhoC) + kD*rhoD)*
           tti^2*(Cos[tt]*(8*xCi + 3*bbCi^4*(1 - 5*bbCi^2*pCi^2)*Cos[tti]) + 
             8*yCi*Sin[tt] + bbCi^4*(7*(-1 + bbCi^2*pCi^2)*Cos[tt - 3*tti] + 
               3*(3 - 7*bbCi^2*pCi^2)*Sin[tt]*Sin[tti]))^2 - 
          2304*bbCi^11*j2^3*(1 - 5*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*
           Sin[tt]*(4*(-(bbCi*dyi) - 8*dbbi*yCi + 60*bbCi^2*dbbi*pCi^2*yCi + 
              5*bbCi^3*pCi*(dyi*pCi + 4*dpi*yCi)) - 3*bbCi^4*
             (bbCi^3*dpi*pCi*(-37 + 105*bbCi^2*pCi^2) + dbbi*(18 - 
                169*bbCi^2*pCi^2 + 315*bbCi^4*pCi^4))*Sin[tti] + 
            7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 15*bbCi^2*pCi^2) + 
              dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*pCi^4))*Sin[3*tti]) - 
          108*bbCi^10*j2*(j2 - 5*bbCi^2*j2*pCi^2)^2*tti*
           (Cos[tt]*(8*xCi + 3*bbCi^4*(1 - 5*bbCi^2*pCi^2)*Cos[tti]) + 
            8*yCi*Sin[tt] + bbCi^4*(7*(-1 + bbCi^2*pCi^2)*Cos[tt - 3*tti] + 
              3*(3 - 7*bbCi^2*pCi^2)*Sin[tt]*Sin[tti]))*
           (32*R^2*(-(kC*rhoC) + kD*rhoD)*(kC*rhoC + kD*rhoD) + 
            bbCi^6*j2^2*(-1 + 5*bbCi^2*pCi^2)*(8*(3*bbCi^5*j2*(-(bbCi*dxi) - 
                  7*dbbi*xCi + 55*bbCi^2*dbbi*pCi^2*xCi + 5*bbCi^3*pCi*
                   (dxi*pCi + 4*dpi*xCi)) + 8*R*(kC*rhoC - kD*rhoD)*yCi)*Cos[
                tt] + 8*(8*R*(-(kC*rhoC) + kD*rhoD)*xCi + 3*bbCi^5*j2*
                 (-(bbCi*dyi) - 7*dbbi*yCi + 55*bbCi^2*dbbi*pCi^2*yCi + 
                  5*bbCi^3*pCi*(dyi*pCi + 4*dpi*yCi)))*Sin[tt] - 
              3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti]*(3*bbCi^5*j2*
                 (30*bbCi^3*dpi*pCi + dbbi*(-11 + 85*bbCi^2*pCi^2))*Cos[tt] + 
                8*R*(-(kC*rhoC) + kD*rhoD)*Sin[tt]) + 7*bbCi^4*Cos[3*tti]*(
                3*bbCi^5*j2*(-22*bbCi^3*dpi*pCi + 30*bbCi^5*dpi*pCi^3 + 
                  dbbi*(11 - 88*bbCi^2*pCi^2 + 85*bbCi^4*pCi^4))*Cos[tt] + 
                8*(-1 + bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[tt]) - 
              3*bbCi^4*(-8*(-3 + 7*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*
                 Cos[tt] + 3*bbCi^5*j2*(-74*bbCi^3*dpi*pCi + 210*bbCi^5*dpi*
                   pCi^3 + dbbi*(33 - 316*bbCi^2*pCi^2 + 595*bbCi^4*pCi^4))*
                 Sin[tt])*Sin[tti] + 7*bbCi^4*(-8*(-1 + bbCi^2*pCi^2)*R*
                 (-(kC*rhoC) + kD*rhoD)*Cos[tt] + 3*bbCi^5*j2*
                 (-22*bbCi^3*dpi*pCi + 30*bbCi^5*dpi*pCi^3 + dbbi*
                   (11 - 88*bbCi^2*pCi^2 + 85*bbCi^4*pCi^4))*Sin[tt])*Sin[
                3*tti])) + 576*bbCi^6*j2*R*(-(kC*rhoC) + kD*rhoD)*Cos[tt]*
           (-3*bbCi^4*(j2 - 5*bbCi^2*j2*pCi^2)^2*(2*bbCi*dbbi + 
              R*(kC*rhoC - kD*rhoD)*tti)*(8*xCi - 3*bbCi^4*(-1 + 
                5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[
                3*tti]) + j2*R*(-(kC*rhoC) + kD*rhoD)*
             (3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*xCi - 3*bbCi^4*
                 (-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Cos[3*tti]) - 4*(-1 + 5*bbCi^2*pCi^2)*(
                8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
                7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti]))) + 
          576*bbCi^6*j2*R*(-(kC*rhoC) + kD*rhoD)*Sin[tt]*
           (-3*bbCi^4*(j2 - 5*bbCi^2*j2*pCi^2)^2*(2*bbCi*dbbi + 
              R*(kC*rhoC - kD*rhoD)*tti)*(8*yCi - 3*bbCi^4*(-3 + 
                7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[
                3*tti]) + j2*R*(-(kC*rhoC) + kD*rhoD)*
             (4*(-1 + 5*bbCi^2*pCi^2)*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
                 Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 
              3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*yCi - 3*bbCi^4*
                 (-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Sin[3*tti]))) + 54*bbCi^16*j2^3*
           (j2 - 5*bbCi^2*j2*pCi^2)^2*(2*bbCi*dbbi + R*(kC*rhoC - kD*rhoD)*
             tti)*Cos[tt]*Sin[tt]*(-((8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*
                Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti])*
              (3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*xCi - 3*bbCi^4*
                  (-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                   bbCi^2*pCi^2)*Cos[3*tti]) - 4*(-1 + 5*bbCi^2*pCi^2)*
                (8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
                 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti]))) - 
            (8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
              7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti])*
             (4*(-1 + 5*bbCi^2*pCi^2)*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
                 Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 
              3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*yCi - 3*bbCi^4*
                 (-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Sin[3*tti]))) - 18*bbCi^12*j2*Cos[tt]*Sin[tt]*
           (-((j2 - 5*bbCi^2*j2*pCi^2)^2*(8*yCi - 3*bbCi^4*(-3 + 
                 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                Sin[3*tti])*(-18*bbCi^9*j2^3*(1 - 5*bbCi^2*pCi^2)*tti*
                (4*(-(bbCi*dxi) - 8*dbbi*xCi + 60*bbCi^2*dbbi*pCi^2*xCi + 
                   5*bbCi^3*pCi*(dxi*pCi + 4*dpi*xCi)) - 9*bbCi^4*
                  (-1 + 5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 
                     15*bbCi^2*pCi^2))*Cos[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*
                    (-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 
                     45*bbCi^4*pCi^4))*Cos[3*tti]) + 8*j2*(12*bbCi^6*dyi*j2 - 
                 60*bbCi^8*dyi*j2*pCi^2 + 8*kC*R*rhoC*xCi - 8*kD*R*rhoD*xCi + 
                 48*bbCi^5*dbbi*j2*yCi - 120*bbCi^8*dpi*j2*pCi*yCi - 
                 360*bbCi^7*dbbi*j2*pCi^2*yCi + 12*bbCi^4*(-1 + 5*bbCi^2*
                    pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Cos[tt] - 28*bbCi^4*
                  (-1 + bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Cos[3*tt] + 
                 15*bbCi^4*kC*R*rhoC*Cos[tti] - 75*bbCi^6*kC*pCi^2*R*rhoC*
                  Cos[tti] - 15*bbCi^4*kD*R*rhoD*Cos[tti] + 75*bbCi^6*kD*
                  pCi^2*R*rhoD*Cos[tti] - 7*bbCi^4*kC*R*rhoC*Cos[3*tti] + 
                 7*bbCi^6*kC*pCi^2*R*rhoC*Cos[3*tti] + 7*bbCi^4*kD*R*rhoD*
                  Cos[3*tti] - 7*bbCi^6*kD*pCi^2*R*rhoD*Cos[3*tti] + 
                 108*bbCi^9*dbbi*j2*Sin[tti] - 198*bbCi^12*dpi*j2*pCi*
                  Sin[tti] - 990*bbCi^11*dbbi*j2*pCi^2*Sin[tti] + 
                 630*bbCi^14*dpi*j2*pCi^3*Sin[tti] + 1890*bbCi^13*dbbi*j2*
                  pCi^4*Sin[tti] - 84*bbCi^9*dbbi*j2*Sin[3*tti] + 
                 126*bbCi^12*dpi*j2*pCi*Sin[3*tti] + 630*bbCi^11*dbbi*j2*
                  pCi^2*Sin[3*tti] - 210*bbCi^14*dpi*j2*pCi^3*Sin[3*tti] - 
                 630*bbCi^13*dbbi*j2*pCi^4*Sin[3*tti]))) + 
            6*bbCi^5*j2^4*(1 - 5*bbCi^2*pCi^2)*(4*(-(bbCi*dxi) - 8*dbbi*xCi + 
                60*bbCi^2*dbbi*pCi^2*xCi + 5*bbCi^3*pCi*(dxi*pCi + 
                  4*dpi*xCi)) - 9*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*
                 pCi + dbbi*(-2 + 15*bbCi^2*pCi^2))*Cos[tti] + 
              7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 15*bbCi^2*pCi^2) + 
                dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*pCi^4))*Cos[3*tti])*
             (4*(-1 + 5*bbCi^2*pCi^2)*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
                 Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 
              3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*yCi - 3*bbCi^4*
                 (-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Sin[3*tti]))) + 36*bbCi^12*j2^2*R*
           (-(kC*rhoC) + kD*rhoD)*Cos[tt]*Sin[tt]*
           (j2^2*(3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*xCi - 
                3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*
                 (-1 + bbCi^2*pCi^2)*Cos[3*tti]) - 4*(-1 + 5*bbCi^2*pCi^2)*(
                8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
                7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti]))*
             (4*(-1 + 5*bbCi^2*pCi^2)*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
                 Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 
              3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*yCi - 3*bbCi^4*
                 (-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Sin[3*tti])) - 64*(j2 - 5*bbCi^2*j2*pCi^2)^2*
             (8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
              7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti])*
             (xCi - (9*bbCi^8*j2^2*(1 - 5*bbCi^2*pCi^2)^2*tti^2*
                (8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
                 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]))/256 + 
              (bbCi^4*(3*(-1 + 5*bbCi^2*pCi^2)*Cos[tt] + 7*Cos[3*tt] + 
                 3*Cos[tti] - 7*Cos[3*tti] + bbCi^2*pCi^2*(-7*Cos[3*tt] - 
                   15*Cos[tti] + 7*Cos[3*tti])))/8 + (3*bbCi^4*j2*
                (-1 + 5*bbCi^2*pCi^2)*tti*(8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*
                    pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                  Sin[3*tti]))/32) - 64*(j2 - 5*bbCi^2*j2*pCi^2)^2*
             (8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
              7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti])*
             (yCi - (3*bbCi^4*j2*(-1 + 5*bbCi^2*pCi^2)*tti*(8*xCi - 
                 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*
                  (-1 + bbCi^2*pCi^2)*Cos[3*tti]))/32 + (bbCi^4*(3*Sin[tt] + 
                 7*(-1 + bbCi^2*pCi^2)*Sin[tt]^3 - 3*Sin[tti] + 
                 (7 - 7*bbCi^2*pCi^2)*Sin[tti]^3))/2 - (9*bbCi^8*j2^2*
                (1 - 5*bbCi^2*pCi^2)^2*tti^2*(8*yCi - 3*bbCi^4*(-3 + 
                   7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                  Sin[3*tti]))/256)) + 9*bbCi^12*j2^2*Sin[tt]^2*
           (6*bbCi^4*j2^3*(1 - 5*bbCi^2*pCi^2)^2*(2*bbCi*dbbi + 
              R*(kC*rhoC - kD*rhoD)*tti)*(8*yCi - 3*bbCi^4*(-3 + 
                7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[
                3*tti])*(-4*(-1 + 5*bbCi^2*pCi^2)*(8*xCi - 3*bbCi^4*
                 (-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Cos[3*tti]) - 3*bbCi^4*j2*(1 - 5*bbCi^2*
                  pCi^2)^2*tti*(8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*
                 Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti])) + 
            R*(-(kC*rhoC) + kD*rhoD)*(32*j2^2*(1 - 5*bbCi^2*pCi^2)^2*
               (8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
                 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti])^2 + 
              128*bbCi^4*(j2 - 5*bbCi^2*j2*pCi^2)^2*(3*Sin[tt] + 
                7*(-1 + bbCi^2*pCi^2)*Sin[tt]^3 - 3*Sin[tti] + 
                (7 - 7*bbCi^2*pCi^2)*Sin[tti]^3)*(-8*yCi + 3*bbCi^4*
                 (-3 + 7*bbCi^2*pCi^2)*Sin[tti] - 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Sin[3*tti]) + 8*(j2 - 5*bbCi^2*j2*pCi^2)^2*(
                8*(-9*bbCi^4*j2*tti*xCi + 45*bbCi^6*j2*pCi^2*tti*xCi - 
                  4*yCi) - 27*bbCi^8*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*Cos[tti] + 
                63*bbCi^8*j2*(1 - 6*bbCi^2*pCi^2 + 5*bbCi^4*pCi^4)*tti*
                 Cos[3*tti])*(8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*
                 Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti]) + 
              27*bbCi^8*(j2 - 5*bbCi^2*j2*pCi^2)^4*tti^2*(8*yCi - 
                 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*
                  (-1 + bbCi^2*pCi^2)*Sin[3*tti])^2)) - 
          18*bbCi^12*j2*Cos[tt]^2*(-((j2 - 5*bbCi^2*j2*pCi^2)^2*
              (8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*
                (-1 + bbCi^2*pCi^2)*Cos[3*tti])*(-18*bbCi^9*j2^3*
                (1 - 5*bbCi^2*pCi^2)*tti*(4*(-(bbCi*dxi) - 8*dbbi*xCi + 
                   60*bbCi^2*dbbi*pCi^2*xCi + 5*bbCi^3*pCi*(dxi*pCi + 
                     4*dpi*xCi)) - 9*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
                  (5*bbCi^3*dpi*pCi + dbbi*(-2 + 15*bbCi^2*pCi^2))*Cos[tti] + 
                 7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 15*bbCi^2*pCi^2) + 
                   dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*pCi^4))*
                  Cos[3*tti]) + 8*j2*(12*bbCi^6*dyi*j2 - 60*bbCi^8*dyi*j2*
                  pCi^2 + 8*kC*R*rhoC*xCi - 8*kD*R*rhoD*xCi + 48*bbCi^5*dbbi*
                  j2*yCi - 120*bbCi^8*dpi*j2*pCi*yCi - 360*bbCi^7*dbbi*j2*
                  pCi^2*yCi + 12*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*R*(-(kC*rhoC) + 
                   kD*rhoD)*Cos[tt] - 28*bbCi^4*(-1 + bbCi^2*pCi^2)*R*
                  (-(kC*rhoC) + kD*rhoD)*Cos[3*tt] + 15*bbCi^4*kC*R*rhoC*
                  Cos[tti] - 75*bbCi^6*kC*pCi^2*R*rhoC*Cos[tti] - 
                 15*bbCi^4*kD*R*rhoD*Cos[tti] + 75*bbCi^6*kD*pCi^2*R*rhoD*
                  Cos[tti] - 7*bbCi^4*kC*R*rhoC*Cos[3*tti] + 7*bbCi^6*kC*
                  pCi^2*R*rhoC*Cos[3*tti] + 7*bbCi^4*kD*R*rhoD*Cos[3*tti] - 
                 7*bbCi^6*kD*pCi^2*R*rhoD*Cos[3*tti] + 108*bbCi^9*dbbi*j2*
                  Sin[tti] - 198*bbCi^12*dpi*j2*pCi*Sin[tti] - 990*bbCi^11*
                  dbbi*j2*pCi^2*Sin[tti] + 630*bbCi^14*dpi*j2*pCi^3*
                  Sin[tti] + 1890*bbCi^13*dbbi*j2*pCi^4*Sin[tti] - 
                 84*bbCi^9*dbbi*j2*Sin[3*tti] + 126*bbCi^12*dpi*j2*pCi*
                  Sin[3*tti] + 630*bbCi^11*dbbi*j2*pCi^2*Sin[3*tti] - 
                 210*bbCi^14*dpi*j2*pCi^3*Sin[3*tti] - 630*bbCi^13*dbbi*j2*
                  pCi^4*Sin[3*tti]))) - 6*bbCi^5*j2^4*(1 - 5*bbCi^2*pCi^2)^2*
             (4*(-(bbCi*dxi) - 8*dbbi*xCi + 60*bbCi^2*dbbi*pCi^2*xCi + 
                5*bbCi^3*pCi*(dxi*pCi + 4*dpi*xCi)) - 9*bbCi^4*(-1 + 
                5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 15*bbCi^2*
                   pCi^2))*Cos[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 
                  15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*
                   pCi^4))*Cos[3*tti])*(-9*bbCi^8*j2*(1 - 5*bbCi^2*pCi^2)^
                2*tti*Cos[tti] + 21*bbCi^8*j2*(1 - 6*bbCi^2*pCi^2 + 
                5*bbCi^4*pCi^4)*tti*Cos[3*tti] - 4*(6*bbCi^4*j2*tti*xCi - 
                30*bbCi^6*j2*pCi^2*tti*xCi + 8*yCi - 3*bbCi^4*(-3 + 
                  7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                 Sin[3*tti]))) - 18*bbCi^12*j2*Cos[tt]*Sin[tt]*
           (-6*bbCi^5*j2^4*(1 - 5*bbCi^2*pCi^2)^2*(4*(-(bbCi*dyi) - 
                8*dbbi*yCi + 60*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*
                 (dyi*pCi + 4*dpi*yCi)) - 3*bbCi^4*(bbCi^3*dpi*pCi*
                 (-37 + 105*bbCi^2*pCi^2) + dbbi*(18 - 169*bbCi^2*pCi^2 + 
                  315*bbCi^4*pCi^4))*Sin[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*
                 (-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 
                  45*bbCi^4*pCi^4))*Sin[3*tti])*(-9*bbCi^8*j2*
               (1 - 5*bbCi^2*pCi^2)^2*tti*Cos[tti] + 21*bbCi^8*j2*(1 - 
                6*bbCi^2*pCi^2 + 5*bbCi^4*pCi^4)*tti*Cos[3*tti] - 
              4*(6*bbCi^4*j2*tti*xCi - 30*bbCi^6*j2*pCi^2*tti*xCi + 8*yCi - 
                3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*
                 (-1 + bbCi^2*pCi^2)*Sin[3*tti])) - (j2 - 5*bbCi^2*j2*pCi^2)^
              2*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
              7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti])*
             (-18*bbCi^9*j2^3*(1 - 5*bbCi^2*pCi^2)*tti*(4*(-(bbCi*dyi) - 
                  8*dbbi*yCi + 60*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*
                   (dyi*pCi + 4*dpi*yCi)) - 3*bbCi^4*(bbCi^3*dpi*pCi*
                   (-37 + 105*bbCi^2*pCi^2) + dbbi*(18 - 169*bbCi^2*pCi^2 + 
                    315*bbCi^4*pCi^4))*Sin[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*
                   (-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 
                    45*bbCi^4*pCi^4))*Sin[3*tti]) - 8*j2*(12*bbCi^6*dxi*j2 - 
                60*bbCi^8*dxi*j2*pCi^2 + 48*bbCi^5*dbbi*j2*xCi - 
                120*bbCi^8*dpi*j2*pCi*xCi - 360*bbCi^7*dbbi*j2*pCi^2*xCi - 
                8*kC*R*rhoC*yCi + 8*kD*R*rhoD*yCi + 18*bbCi^9*j2*
                 (-1 + 5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 
                    15*bbCi^2*pCi^2))*Cos[tti] - 42*bbCi^9*j2*
                 (bbCi^3*dpi*pCi*(-3 + 5*bbCi^2*pCi^2) + dbbi*(2 - 15*bbCi^2*
                     pCi^2 + 15*bbCi^4*pCi^4))*Cos[3*tti] - 36*bbCi^4*kC*R*
                 rhoC*Sin[tt] + 84*bbCi^6*kC*pCi^2*R*rhoC*Sin[tt] + 
                36*bbCi^4*kD*R*rhoD*Sin[tt] - 84*bbCi^6*kD*pCi^2*R*rhoD*
                 Sin[tt] + 28*bbCi^4*kC*R*rhoC*Sin[3*tt] - 28*bbCi^6*kC*pCi^2*
                 R*rhoC*Sin[3*tt] - 28*bbCi^4*kD*R*rhoD*Sin[3*tt] + 
                28*bbCi^6*kD*pCi^2*R*rhoD*Sin[3*tt] - 21*bbCi^4*kC*R*rhoC*
                 Sin[tti] + 81*bbCi^6*kC*pCi^2*R*rhoC*Sin[tti] + 
                21*bbCi^4*kD*R*rhoD*Sin[tti] - 81*bbCi^6*kD*pCi^2*R*rhoD*
                 Sin[tti] + 7*bbCi^4*kC*R*rhoC*Sin[3*tti] - 7*bbCi^6*kC*pCi^2*
                 R*rhoC*Sin[3*tti] - 7*bbCi^4*kD*R*rhoD*Sin[3*tti] + 
                7*bbCi^6*kD*pCi^2*R*rhoD*Sin[3*tti]))) - 
          18*bbCi^12*j2*Sin[tt]^2*(6*bbCi^5*j2^4*(1 - 5*bbCi^2*pCi^2)*
             (4*(-(bbCi*dyi) - 8*dbbi*yCi + 60*bbCi^2*dbbi*pCi^2*yCi + 
                5*bbCi^3*pCi*(dyi*pCi + 4*dpi*yCi)) - 3*bbCi^4*(bbCi^3*dpi*
                 pCi*(-37 + 105*bbCi^2*pCi^2) + dbbi*(18 - 169*bbCi^2*pCi^2 + 
                  315*bbCi^4*pCi^4))*Sin[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*
                 (-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 
                  45*bbCi^4*pCi^4))*Sin[3*tti])*(4*(-1 + 5*bbCi^2*pCi^2)*(
                8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
                7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 3*bbCi^4*j2*
               (1 - 5*bbCi^2*pCi^2)^2*tti*(8*yCi - 3*bbCi^4*(-3 + 
                  7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                 Sin[3*tti])) - (j2 - 5*bbCi^2*j2*pCi^2)^2*
             (8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
              7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti])*
             (-18*bbCi^9*j2^3*(1 - 5*bbCi^2*pCi^2)*tti*(4*(-(bbCi*dyi) - 
                  8*dbbi*yCi + 60*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*
                   (dyi*pCi + 4*dpi*yCi)) - 3*bbCi^4*(bbCi^3*dpi*pCi*
                   (-37 + 105*bbCi^2*pCi^2) + dbbi*(18 - 169*bbCi^2*pCi^2 + 
                    315*bbCi^4*pCi^4))*Sin[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*
                   (-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 
                    45*bbCi^4*pCi^4))*Sin[3*tti]) - 8*j2*(12*bbCi^6*dxi*j2 - 
                60*bbCi^8*dxi*j2*pCi^2 + 48*bbCi^5*dbbi*j2*xCi - 
                120*bbCi^8*dpi*j2*pCi*xCi - 360*bbCi^7*dbbi*j2*pCi^2*xCi - 
                8*kC*R*rhoC*yCi + 8*kD*R*rhoD*yCi + 18*bbCi^9*j2*
                 (-1 + 5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 
                    15*bbCi^2*pCi^2))*Cos[tti] - 42*bbCi^9*j2*
                 (bbCi^3*dpi*pCi*(-3 + 5*bbCi^2*pCi^2) + dbbi*(2 - 15*bbCi^2*
                     pCi^2 + 15*bbCi^4*pCi^4))*Cos[3*tti] - 36*bbCi^4*kC*R*
                 rhoC*Sin[tt] + 84*bbCi^6*kC*pCi^2*R*rhoC*Sin[tt] + 
                36*bbCi^4*kD*R*rhoD*Sin[tt] - 84*bbCi^6*kD*pCi^2*R*rhoD*
                 Sin[tt] + 28*bbCi^4*kC*R*rhoC*Sin[3*tt] - 28*bbCi^6*kC*pCi^2*
                 R*rhoC*Sin[3*tt] - 28*bbCi^4*kD*R*rhoD*Sin[3*tt] + 
                28*bbCi^6*kD*pCi^2*R*rhoD*Sin[3*tt] - 21*bbCi^4*kC*R*rhoC*
                 Sin[tti] + 81*bbCi^6*kC*pCi^2*R*rhoC*Sin[tti] + 
                21*bbCi^4*kD*R*rhoD*Sin[tti] - 81*bbCi^6*kD*pCi^2*R*rhoD*
                 Sin[tti] + 7*bbCi^4*kC*R*rhoC*Sin[3*tti] - 7*bbCi^6*kC*pCi^2*
                 R*rhoC*Sin[3*tti] - 7*bbCi^4*kD*R*rhoD*Sin[3*tti] + 
                7*bbCi^6*kD*pCi^2*R*rhoD*Sin[3*tti]))) + 
          9*bbCi^12*j2^2*Cos[tt]^2*(6*bbCi^4*j2^3*(1 - 5*bbCi^2*pCi^2)^2*
             (2*bbCi*dbbi + R*(kC*rhoC - kD*rhoD)*tti)*(8*xCi - 
              3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                bbCi^2*pCi^2)*Cos[3*tti])*(-3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^
                2*tti*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
                7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 
              4*(-1 + 5*bbCi^2*pCi^2)*(8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*
                 Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti])) + 
            R*(-(kC*rhoC) + kD*rhoD)*(-256*(j2 - 5*bbCi^2*j2*pCi^2)^2*xCi*(
                8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
                7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 
              27*bbCi^8*(j2 - 5*bbCi^2*j2*pCi^2)^4*tti^2*(8*xCi - 
                 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*
                  (-1 + bbCi^2*pCi^2)*Cos[3*tti])^2 + 32*j2^2*
               (1 - 5*bbCi^2*pCi^2)^2*(8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*
                  Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti])^2 - 
              8*bbCi^4*(j2 - 5*bbCi^2*j2*pCi^2)^2*(8*xCi - 3*bbCi^4*
                 (-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                  bbCi^2*pCi^2)*Cos[3*tti])*(-72*j2*tti*yCi + 360*bbCi^2*j2*
                 pCi^2*tti*yCi + 12*(-1 + 5*bbCi^2*pCi^2)*Cos[tt] - 
                28*(-1 + bbCi^2*pCi^2)*Cos[3*tt] + 12*Cos[tti] - 
                60*bbCi^2*pCi^2*Cos[tti] - 28*Cos[3*tti] + 28*bbCi^2*pCi^2*
                 Cos[3*tti] - 81*bbCi^4*j2*tti*Sin[tti] + 594*bbCi^6*j2*pCi^2*
                 tti*Sin[tti] - 945*bbCi^8*j2*pCi^4*tti*Sin[tti] + 
                63*bbCi^4*j2*tti*Sin[3*tti] - 378*bbCi^6*j2*pCi^2*tti*
                 Sin[3*tti] + 315*bbCi^8*j2*pCi^4*tti*Sin[3*tti]))) + 
          32*bbCi^2*j2*R*(72*bbCi^9*j2^2*kC*(1 - 5*bbCi^2*pCi^2)*rhoC*Cos[tt]*
             (4*(-(bbCi*dxi) - 8*dbbi*xCi + 60*bbCi^2*dbbi*pCi^2*xCi + 
                5*bbCi^3*pCi*(dxi*pCi + 4*dpi*xCi)) - 9*bbCi^4*(-1 + 
                5*bbCi^2*pCi^2)*(5*bbCi^3*dpi*pCi + dbbi*(-2 + 15*bbCi^2*
                   pCi^2))*Cos[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 
                  15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*
                   pCi^4))*Cos[3*tti]) + 72*bbCi^9*j2^2*kC*
             (1 - 5*bbCi^2*pCi^2)*rhoC*Sin[tt]*(4*(-(bbCi*dyi) - 8*dbbi*yCi + 
                60*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*(dyi*pCi + 
                  4*dpi*yCi)) - 3*bbCi^4*(bbCi^3*dpi*pCi*(-37 + 105*bbCi^2*
                   pCi^2) + dbbi*(18 - 169*bbCi^2*pCi^2 + 315*bbCi^4*pCi^4))*
               Sin[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 15*bbCi^2*pCi^2) + 
                dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*pCi^4))*Sin[3*tti]) + 
            3*bbCi^4*(-(kC*rhoC) + kD*rhoD)*Sin[tt]*
             (3*bbCi^4*(j2 - 5*bbCi^2*j2*pCi^2)^2*(-4*bbCi*dbbi + 
                3*bbCi^6*j2 - 15*bbCi^8*j2*pCi^2 - 4*kC*R*rhoC*tti + 
                18*bbCi^6*j2*(-1 + bbCi^2*pCi^2)*Cos[2*tt])*(8*yCi - 
                3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*
                 (-1 + bbCi^2*pCi^2)*Sin[3*tti]) + 4*j2*kC*(1 - 
                5*bbCi^2*pCi^2)*R*rhoC*(4*(8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*
                     pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                   Cos[3*tti]) - 3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)*tti*
                 (8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
                  7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti]))) + 
            Cos[tt]*(18*bbCi^8*kC*(j2 - 5*bbCi^2*j2*pCi^2)^2*rhoC*(
                2*bbCi*dbbi + R*(kC*rhoC - kD*rhoD)*tti)*(8*xCi - 
                3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*
                 (-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 9*bbCi^9*(j2 - 
                 5*bbCi^2*j2*pCi^2)^2*(-3*bbCi^5*j2*kC*rhoC + 15*bbCi^7*j2*kC*
                 pCi^2*rhoC - 4*dbbi*kD*rhoD + 3*bbCi^5*j2*kD*rhoD - 
                15*bbCi^7*j2*kD*pCi^2*rhoD + 15*bbCi^5*j2*(-1 + bbCi^2*pCi^2)*
                 (-(kC*rhoC) + kD*rhoD)*Cos[2*tt] + 3*bbCi^5*j2*
                 (-1 + bbCi^2*pCi^2)*(-(kC*rhoC) + kD*rhoD)*Cos[2*tti])*(
                8*xCi - 3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 
                7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[3*tti]) - (-(kC*rhoC) + 
                kD*rhoD)*(9*bbCi^8*(j2 - 5*bbCi^2*j2*pCi^2)^2*
                 (2*kC*R*rhoC*tti - 3*bbCi^6*j2*(-1 + bbCi^2*pCi^2)*
                   (Cos[2*tt] - Cos[2*tti]))*(8*xCi - 3*bbCi^4*(-1 + 
                    5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                   Cos[3*tti]) + 12*bbCi^4*j2*kC*(1 - 5*bbCi^2*pCi^2)*R*rhoC*
                 (3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)*tti*(8*xCi - 3*bbCi^4*
                     (-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                      bbCi^2*pCi^2)*Cos[3*tti]) + 4*(8*yCi - 3*bbCi^4*
                     (-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + 
                      bbCi^2*pCi^2)*Sin[3*tti]))))) + 
          128*R*(-9*bbCi^10*j2^3*(1 - 5*bbCi^2*pCi^2)^2*(-(kC*rhoC) + 
              kD*rhoD)*(-2*kC*R*rhoC*tti + 3*bbCi^6*j2*(-1 + bbCi^2*pCi^2)*
               Cos[2*tt] - 3*bbCi^6*j2*(-1 + bbCi^2*pCi^2)*Cos[2*tti])*
             (7*bbCi^4*(-1 + bbCi^2*pCi^2)*Cos[tt - 3*tti] + Cos[tt]*(8*xCi - 
                3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti]) + Sin[tt]*(8*yCi + 
                9*bbCi^4*Sin[tti] - 21*bbCi^6*pCi^2*Sin[tti])) + 
            4*kC*rhoC*(96*R^2*(-(kC*rhoC) + kD*rhoD)^2 - 9*bbCi^11*j2^3*(1 - 
                5*bbCi^2*pCi^2)*Cos[tt]*(4*(-(bbCi*dxi) - 8*dbbi*xCi + 
                  60*bbCi^2*dbbi*pCi^2*xCi + 5*bbCi^3*pCi*(dxi*pCi + 
                    4*dpi*xCi)) - 9*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*
                 (5*bbCi^3*dpi*pCi + dbbi*(-2 + 15*bbCi^2*pCi^2))*Cos[tti] + 
                7*bbCi^4*(bbCi^3*dpi*pCi*(-11 + 15*bbCi^2*pCi^2) + 
                  dbbi*(6 - 47*bbCi^2*pCi^2 + 45*bbCi^4*pCi^4))*Cos[3*tti]) - 
              9*bbCi^11*j2^3*(1 - 5*bbCi^2*pCi^2)*Sin[tt]*(4*(-(bbCi*dyi) - 
                  8*dbbi*yCi + 60*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*
                   (dyi*pCi + 4*dpi*yCi)) - 3*bbCi^4*(bbCi^3*dpi*pCi*
                   (-37 + 105*bbCi^2*pCi^2) + dbbi*(18 - 169*bbCi^2*pCi^2 + 
                    315*bbCi^4*pCi^4))*Sin[tti] + 7*bbCi^4*(bbCi^3*dpi*pCi*
                   (-11 + 15*bbCi^2*pCi^2) + dbbi*(6 - 47*bbCi^2*pCi^2 + 
                    45*bbCi^4*pCi^4))*Sin[3*tti]) + (bbCi^2*j2*Cos[tt]*
                (-9*bbCi^8*(j2 - 5*bbCi^2*j2*pCi^2)^2*(2*bbCi*dbbi + 
                   R*(kC*rhoC - kD*rhoD)*tti)*(8*xCi - 3*bbCi^4*(-1 + 
                     5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + bbCi^2*pCi^2)*
                    Cos[3*tti]) + 6*bbCi^4*j2*R*(-(kC*rhoC) + kD*rhoD)*
                  (3*bbCi^4*j2*(1 - 5*bbCi^2*pCi^2)^2*tti*(8*xCi - 3*bbCi^4*
                      (-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*(-1 + 
                       bbCi^2*pCi^2)*Cos[3*tti]) - 4*(-1 + 5*bbCi^2*pCi^2)*
                    (8*yCi - 3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 
                     7*bbCi^4*(-1 + bbCi^2*pCi^2)*Sin[3*tti]))))/2 + 
              (bbCi^2*j2*Sin[tt]*(-9*bbCi^8*(j2 - 5*bbCi^2*j2*pCi^2)^2*
                  (2*bbCi*dbbi + R*(kC*rhoC - kD*rhoD)*tti)*(8*yCi - 
                   3*bbCi^4*(-3 + 7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*
                    (-1 + bbCi^2*pCi^2)*Sin[3*tti]) + 6*bbCi^4*j2*R*
                  (-(kC*rhoC) + kD*rhoD)*(4*(-1 + 5*bbCi^2*pCi^2)*(8*xCi - 
                     3*bbCi^4*(-1 + 5*bbCi^2*pCi^2)*Cos[tti] + 7*bbCi^4*
                      (-1 + bbCi^2*pCi^2)*Cos[3*tti]) + 3*bbCi^4*j2*
                    (1 - 5*bbCi^2*pCi^2)^2*tti*(8*yCi - 3*bbCi^4*(-3 + 
                       7*bbCi^2*pCi^2)*Sin[tti] + 7*bbCi^4*(-1 + bbCi^2*
                        pCi^2)*Sin[3*tti]))))/2))))/bbCi - 
       1024*bbCi^2*(6*(-96*bbCi*dbbi^2 + 128*dbbi^3 - 48*bbCi^6*dbbi*j2 - 
           24*bbCi^5*dbbi^2*j2 - 159*bbCi^10*dbbi*j2^2 + 118*bbCi^13*dpi*j2^2*
            pCi + 432*bbCi^7*dbbi^2*j2*pCi^2 + 472*bbCi^12*dbbi*j2^2*pCi^2 + 
           14*bbCi^15*dpi*j2^2*pCi^3 + 35*bbCi^14*dbbi*j2^2*pCi^4 + 
           288*bbCi^8*dbbi*j2*pCi*(2*dpi + pCi) + 72*bbCi^9*dpi*j2*
            (dpi + 2*pCi) - 32*bbCi^3*j2^2*(dxi*xCi + dyi*yCi) + 
           32*bbCi^2*dbbi*(2 + j2^2*(xCi^2 + yCi^2))) + 
         9*bbCi^10*j2^2*(2*bbCi^3*dpi*pCi*(-48 + 55*bbCi^2*pCi^2) + 
           dbbi*(123 - 384*bbCi^2*pCi^2 + 275*bbCi^4*pCi^4))*Cos[tt]*
          Cos[5*tti] - 147*bbCi^10*j2^2*(-1 + bbCi^2*pCi^2)*
          (2*bbCi^3*dpi*pCi + dbbi*(-3 + 5*bbCi^2*pCi^2))*Cos[2*tt]*
          Cos[6*tti] - 6*bbCi^5*j2*Cos[4*tti]*
          (51*bbCi^5*j2*(-1 + bbCi^2*pCi^2)*(2*bbCi^3*dpi*pCi + 
             dbbi*(-3 + 5*bbCi^2*pCi^2)) - 18*bbCi*j2*(-(bbCi*dxi) - 
             2*dbbi*xCi + 4*bbCi^2*dbbi*pCi^2*xCi + bbCi^3*pCi*
              (dxi*pCi + 2*dpi*xCi))*Cos[tt] - 42*bbCi^5*j2*
            (-4*bbCi^3*dpi*pCi + 6*bbCi^5*dpi*pCi^3 + 
             dbbi*(3 - 16*bbCi^2*pCi^2 + 15*bbCi^4*pCi^4))*Cos[2*tt] + 
           18*bbCi*j2*(-(bbCi*dyi) - 2*dbbi*yCi + 4*bbCi^2*dbbi*pCi^2*yCi + 
             bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi))*Sin[tt] - 
           28*(-1 + bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[2*tt]) - 
         3*bbCi*j2*Cos[2*tti]*(24*bbCi^4*(-8*bbCi*dbbi - 4*dbbi^2 + 
             4*bbCi^4*dpi^2 - 3*bbCi^5*dbbi*j2 + 2*bbCi^3*dpi*
              (4*bbCi + 16*dbbi - bbCi^5*j2)*pCi + 8*bbCi^2*dbbi*
              (2*bbCi + 3*dbbi - bbCi^5*j2)*pCi^2 + 6*bbCi^10*dpi*j2*pCi^3 + 
             15*bbCi^9*dbbi*j2*pCi^4) - 
           16*(3*bbCi^5*j2*(bbCi*dxi*(-4 + 3*bbCi^2*pCi^2) - 8*dbbi*xCi + 6*
                bbCi^2*pCi*(bbCi*dpi + 2*dbbi*pCi)*xCi) + 
             R*(-(kC*rhoC) + kD*rhoD)*yCi)*Cos[tt] + 3*bbCi^9*j2*
            (2*bbCi^3*dpi*pCi*(-77 + 157*bbCi^2*pCi^2) + 
             dbbi*(135 - 616*bbCi^2*pCi^2 + 785*bbCi^4*pCi^4))*Cos[2*tt] - 
           16*(R*(-(kC*rhoC) + kD*rhoD)*xCi + 3*bbCi^5*j2*(-5*bbCi*dyi - 10*
                dbbi*yCi + 24*bbCi^2*dbbi*pCi^2*yCi + 6*bbCi^3*pCi*
                (dyi*pCi + 2*dpi*yCi)))*Sin[tt] + 48*bbCi^4*
            (-1 + 3*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[2*tt]) - 
         bbCi^5*j2*Cos[3*tti]*(264*bbCi*j2*(-(bbCi*dxi) - 2*dbbi*xCi + 
             4*bbCi^2*dbbi*pCi^2*xCi + bbCi^3*pCi*(dxi*pCi + 2*dpi*xCi)) + 
           12*(28*bbCi*dbbi + 14*dbbi^2 - 14*bbCi^4*dpi^2 + 
             45*bbCi^5*dbbi*j2 - bbCi^3*dpi*(28*bbCi + 112*dbbi + 101*bbCi^5*
                j2)*pCi - 4*bbCi^2*dbbi*(14*bbCi + 21*dbbi + 101*bbCi^5*j2)*
              pCi^2 + 208*bbCi^10*dpi*j2*pCi^3 + 520*bbCi^9*dbbi*j2*pCi^4)*
            Cos[tt] + 168*bbCi*j2*(-(bbCi*dxi) - 2*dbbi*xCi + 
             4*bbCi^2*dbbi*pCi^2*xCi + bbCi^3*pCi*(dxi*pCi + 2*dpi*xCi))*
            Cos[2*tt] - 21*bbCi^5*j2*(2*bbCi^3*dpi*pCi*(-4 + 7*bbCi^2*
                pCi^2) + dbbi*(3 - 32*bbCi^2*pCi^2 + 35*bbCi^4*pCi^4))*
            Cos[3*tt] + 8*(-44 + 53*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*
            Sin[tt] + 168*bbCi*j2*(-(bbCi*dyi) - 2*dbbi*yCi + 
             4*bbCi^2*dbbi*pCi^2*yCi + bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi))*
            Sin[2*tt]) - 3*Cos[tti]*(-24*bbCi^6*j2^2*(3*bbCi*dxi + 
             6*dbbi*xCi + 4*bbCi^2*dbbi*pCi^2*xCi + bbCi^3*pCi*
              (dxi*pCi + 2*dpi*xCi)) + 3*bbCi^5*j2*(-16*bbCi*dbbi - 
             8*dbbi^2 + 40*bbCi^4*dpi^2 - 189*bbCi^5*dbbi*j2 + 
             16*bbCi^3*dpi*(5*bbCi + 20*dbbi + 18*bbCi^5*j2)*pCi + 
             16*bbCi^2*dbbi*(10*bbCi + 15*dbbi + 72*bbCi^5*j2)*pCi^2 - 
             546*bbCi^10*dpi*j2*pCi^3 - 1365*bbCi^9*dbbi*j2*pCi^4)*Cos[tt] - 
           8*bbCi*j2*(3*bbCi^5*j2*(-(bbCi*dxi) - 2*dbbi*xCi + 20*bbCi^2*dbbi*
                pCi^2*xCi + 5*bbCi^3*pCi*(dxi*pCi + 2*dpi*xCi)) + 
             8*R*(kC*rhoC - kD*rhoD)*yCi)*Cos[2*tt] + 3*bbCi^10*j2^2*
            (2*bbCi^3*dpi*pCi*(-6 + 35*bbCi^2*pCi^2) + 
             dbbi*(3 - 48*bbCi^2*pCi^2 + 175*bbCi^4*pCi^4))*Cos[3*tt] + 
           4*R*(-64*dbbi*kD*rhoD + 37*bbCi^5*j2*(kC*rhoC - kD*rhoD) + 
             16*bbCi*(-(kC*rhoC) + kD*rhoD) + 49*bbCi^7*j2*pCi^2*
              (-(kC*rhoC) + kD*rhoD))*Sin[tt] - 8*bbCi*j2*
            (8*R*(-(kC*rhoC) + kD*rhoD)*xCi + 3*bbCi^5*j2*(-(bbCi*dyi) - 2*
                dbbi*yCi + 20*bbCi^2*dbbi*pCi^2*yCi + 5*bbCi^3*pCi*
                (dyi*pCi + 2*dpi*yCi)))*Sin[2*tt] + 4*bbCi^5*j2*
            (-1 + 7*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Sin[3*tt]) - 
         4*bbCi*j2*(6*(-8*bbCi^2*dxi - 5*bbCi^6*dxi*j2 - 10*bbCi^5*dbbi*j2*
              xCi + 104*bbCi^7*dbbi*j2*pCi^2*xCi + 16*bbCi*dbbi*(dxi + xCi) + 
             26*bbCi^8*j2*pCi*(dxi*pCi + 2*dpi*xCi) - 
             6*(4*dbbi^2*xCi + R*(-(kC*rhoC) + kD*rhoD)*yCi))*Cos[tt] + 
           3*bbCi*(8*bbCi^4*dbbi + 4*bbCi^3*dbbi^2 - 15*bbCi^8*dbbi*j2 - 
             24*bbCi^11*dpi*j2*pCi - 24*bbCi^5*dbbi^2*pCi^2 - 
             96*bbCi^10*dbbi*j2*pCi^2 + 58*bbCi^13*dpi*j2*pCi^3 + 
             145*bbCi^12*dbbi*j2*pCi^4 - 16*bbCi^6*dbbi*pCi*(2*dpi + pCi) - 
             4*bbCi^7*dpi*(dpi + 2*pCi) + 16*bbCi*j2*(dxi*xCi - dyi*yCi) + 
             16*dbbi*j2*(-xCi^2 + yCi^2))*Cos[2*tt] + 3*bbCi^5*j2*
            (2*dbbi*xCi + bbCi*(dxi - 7*bbCi^2*dxi*pCi^2 - 14*bbCi*pCi*
                (bbCi*dpi + 2*dbbi*pCi)*xCi))*Cos[3*tt] - 
           6*bbCi^9*j2*(bbCi^3*dpi*pCi*(-5 + 8*bbCi^2*pCi^2) + 
             dbbi*(3 - 20*bbCi^2*pCi^2 + 20*bbCi^4*pCi^4))*Cos[4*tt] + 
           6*(-8*bbCi^2*dyi - 13*bbCi^6*dyi*j2 - 26*bbCi^5*dbbi*j2*yCi + 
             112*bbCi^7*dbbi*j2*pCi^2*yCi + 16*bbCi*dbbi*(dyi + yCi) - 
             6*(kC*R*rhoC*xCi - kD*R*rhoD*xCi + 4*dbbi^2*yCi) + 
             28*bbCi^8*j2*pCi*(dyi*pCi + 2*dpi*yCi))*Sin[tt] + 
           4*bbCi*(14*bbCi^5*pCi^2*R*(kC*rhoC - kD*rhoD) + 
             5*bbCi^3*R*(-(kC*rhoC) + kD*rhoD) - 24*dbbi*j2*xCi*yCi + 
             12*bbCi*j2*(dyi*xCi + dxi*yCi))*Sin[2*tt] + 
           3*bbCi^5*j2*(2*dbbi*yCi + bbCi*(dyi - 7*bbCi^2*dyi*pCi^2 - 14*bbCi*
                pCi*(bbCi*dpi + 2*dbbi*pCi)*yCi))*Sin[3*tt]) - 
         3*(4*R*(64*dbbi*kD*rhoD + 16*bbCi*(kC*rhoC - kD*rhoD) + 
             25*bbCi^5*j2*(kC*rhoC - kD*rhoD) + 13*bbCi^7*j2*pCi^2*
              (-(kC*rhoC) + kD*rhoD))*Cos[tt] + bbCi*j2*
            (-24*bbCi^5*j2*(-7*bbCi*dyi - 14*dbbi*yCi + 44*bbCi^2*dbbi*pCi^2*
                yCi + 11*bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi)) + 
             8*(8*R*(-(kC*rhoC) + kD*rhoD)*xCi + 3*bbCi^5*j2*(-3*bbCi*dyi - 
                 6*dbbi*yCi + 28*bbCi^2*dbbi*pCi^2*yCi + 7*bbCi^3*pCi*
                  (dyi*pCi + 2*dpi*yCi)))*Cos[2*tt] - 4*bbCi^4*
              (-1 + 7*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*Cos[3*tt] + 
             3*bbCi^4*(-48*bbCi*dbbi - 24*dbbi^2 + 56*bbCi^4*dpi^2 - 213*
                bbCi^5*dbbi*j2 + 16*bbCi^3*dpi*(7*bbCi + 28*dbbi + 
                 5*bbCi^5*j2)*pCi + 16*bbCi^2*dbbi*(14*bbCi + 21*dbbi + 
                 20*bbCi^5*j2)*pCi^2 + 78*bbCi^10*dpi*j2*pCi^3 + 195*bbCi^9*
                dbbi*j2*pCi^4)*Sin[tt] - 8*(3*bbCi^5*j2*(-3*bbCi*dxi - 
                 6*dbbi*xCi + 28*bbCi^2*dbbi*pCi^2*xCi + 7*bbCi^3*pCi*
                  (dxi*pCi + 2*dpi*xCi)) + 8*R*(kC*rhoC - kD*rhoD)*yCi)*
              Sin[2*tt] + 3*bbCi^9*j2*(14*bbCi^3*dpi*pCi*(-2 + 7*bbCi^2*
                  pCi^2) + dbbi*(9 - 112*bbCi^2*pCi^2 + 245*bbCi^4*pCi^4))*
              Sin[3*tt]))*Sin[tti] - 3*bbCi*j2*(8*bbCi^4*(5 + bbCi^2*pCi^2)*R*
            (-(kC*rhoC) + kD*rhoD) + 16*(R*(-(kC*rhoC) + kD*rhoD)*xCi + 
             6*bbCi^5*j2*(-(bbCi*dyi) - 2*dbbi*yCi + 8*bbCi^2*dbbi*pCi^2*
                yCi + 2*bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi)))*Cos[tt] - 
           48*bbCi^4*(-1 + 3*bbCi^2*pCi^2)*R*(-(kC*rhoC) + kD*rhoD)*
            Cos[2*tt] - 16*(3*bbCi^5*j2*(-(bbCi*dxi) - 2*dbbi*xCi + 4*bbCi^2*
                dbbi*pCi^2*xCi + bbCi^3*pCi*(dxi*pCi + 2*dpi*xCi)) + 
             R*(-(kC*rhoC) + kD*rhoD)*yCi)*Sin[tt] + 3*bbCi^9*j2*
            (-94*bbCi^3*dpi*pCi + 238*bbCi^5*dpi*pCi^3 + 
             dbbi*(69 - 376*bbCi^2*pCi^2 + 595*bbCi^4*pCi^4))*Sin[2*tt])*
          Sin[2*tti] + bbCi^5*j2*(8*(-44 + 53*bbCi^2*pCi^2)*R*
            (-(kC*rhoC) + kD*rhoD)*Cos[tt] + 
           3*(-56*bbCi*j2*(bbCi*dyi + 2*dbbi*yCi - 4*bbCi^2*dbbi*pCi^2*yCi - 
               bbCi^3*pCi*(dyi*pCi + 2*dpi*yCi))*Cos[2*tt] - 
             4*(28*bbCi*dbbi + 2*(7*dbbi^2 - 7*bbCi^4*dpi^2 + 60*bbCi^5*dbbi*
                  j2) - bbCi^3*dpi*(28*bbCi + 112*dbbi + 121*bbCi^5*j2)*pCi - 
               4*bbCi^2*dbbi*(14*bbCi + 21*dbbi + 121*bbCi^5*j2)*pCi^2 + 198*
                bbCi^10*dpi*j2*pCi^3 + 495*bbCi^9*dbbi*j2*pCi^4)*Sin[tt] + 
             bbCi*j2*(176*dbbi*yCi + 88*bbCi*(dyi - bbCi^2*dyi*pCi^2 - 
                 2*bbCi*pCi*(bbCi*dpi + 2*dbbi*pCi)*yCi) + 56*(bbCi*dxi + 
                 2*dbbi*xCi - 4*bbCi^2*dbbi*pCi^2*xCi - bbCi^3*pCi*
                  (dxi*pCi + 2*dpi*xCi))*Sin[2*tt] + 7*bbCi^4*
                (2*bbCi^3*dpi*pCi*(-4 + 7*bbCi^2*pCi^2) + dbbi*
                  (3 - 32*bbCi^2*pCi^2 + 35*bbCi^4*pCi^4))*Sin[3*tt])))*
          Sin[3*tti] + 12*bbCi^5*j2*(-14*(-1 + bbCi^2*pCi^2)*R*
            (-(kC*rhoC) + kD*rhoD)*Cos[2*tt] + 9*bbCi*j2*(-(bbCi*dxi) - 
             2*dbbi*xCi + 4*bbCi^2*dbbi*pCi^2*xCi + bbCi^3*pCi*
              (dxi*pCi + 2*dpi*xCi))*Sin[tt] + 3*bbCi*j2*Cos[tt]*
            (3*(-(bbCi*dyi) - 2*dbbi*yCi + 4*bbCi^2*dbbi*pCi^2*yCi + bbCi^3*
                pCi*(dyi*pCi + 2*dpi*yCi)) + 14*bbCi^4*(2*bbCi^3*dpi*pCi*
                (-2 + 3*bbCi^2*pCi^2) + dbbi*(3 - 16*bbCi^2*pCi^2 + 
                 15*bbCi^4*pCi^4))*Sin[tt]))*Sin[4*tti] + 
         9*bbCi^10*j2^2*(2*bbCi^3*dpi*pCi*(-48 + 55*bbCi^2*pCi^2) + 
           dbbi*(123 - 384*bbCi^2*pCi^2 + 275*bbCi^4*pCi^4))*Sin[tt]*
          Sin[5*tti] - 147*bbCi^10*j2^2*(-1 + bbCi^2*pCi^2)*
          (2*bbCi^3*dpi*pCi + dbbi*(-3 + 5*bbCi^2*pCi^2))*Sin[2*tt]*
          Sin[6*tti]) - 768*bbCi*(tt - tti)*(-256*bbCi^2*kC*R*rhoC + 
         536*bbCi^10*j2^2*kC*R*rhoC + 384*bbCi^8*j2*kC*pCi^2*R*rhoC - 
         2320*bbCi^12*j2^2*kC*pCi^2*R*rhoC + 2936*bbCi^14*j2^2*kC*pCi^4*R*
          rhoC + 256*bbCi^2*kD*R*rhoD - 1024*bbCi*dbbi*kD*R*rhoD + 
         2560*dbbi^2*kD*R*rhoD - 536*bbCi^10*j2^2*kD*R*rhoD - 
         768*bbCi^8*dpi*j2*kD*pCi*R*rhoD - 384*bbCi^8*j2*kD*pCi^2*R*rhoD - 
         768*bbCi^7*dbbi*j2*kD*pCi^2*R*rhoD + 2320*bbCi^12*j2^2*kD*pCi^2*R*
          rhoD - 2936*bbCi^14*j2^2*kD*pCi^4*R*rhoD - 256*bbCi^2*j2^2*kC*R*
          rhoC*xCi^2 + 256*bbCi^2*j2^2*kD*R*rhoD*xCi^2 - 
         256*bbCi^2*j2^2*kC*R*rhoC*yCi^2 + 256*bbCi^2*j2^2*kD*R*rhoD*yCi^2 + 
         4*bbCi*j2*(141*bbCi^11*dyi*j2^2 + 384*dbbi*kD*R*rhoD*xCi - 
           96*bbCi*R*(dxi*kD*rhoD - kC*rhoC*xCi + kD*rhoD*xCi) + 
           846*bbCi^10*dbbi*j2^2*yCi - 8256*bbCi^12*dbbi*j2^2*pCi^2*yCi + 
           13950*bbCi^14*dbbi*j2^2*pCi^4*yCi + 96*bbCi^6*dbbi*j2*
            (dyi + yCi) + 8*bbCi^5*j2*(-7*kC*R*rhoC*xCi + 7*kD*R*rhoD*xCi + 
             6*dbbi^2*yCi) - 1032*bbCi^13*j2^2*pCi*(dyi*pCi + 2*dpi*yCi) + 
           1395*bbCi^15*j2^2*pCi^3*(dyi*pCi + 4*dpi*yCi) - 
           960*bbCi^8*dbbi*j2*pCi*(dyi*pCi + (2*dpi + pCi)*yCi) - 
           240*bbCi^9*j2*(dyi*pCi^2 + dpi^2*yCi + 2*dpi*pCi*(dyi + yCi)) + 
           16*bbCi^7*j2*(3*dyi + pCi^2*(-(kC*R*rhoC*xCi) + kD*R*rhoD*xCi - 90*
                dbbi^2*yCi)))*Cos[tt] + 16*bbCi^2*j2*
          (-186*bbCi^10*j2*pCi^2*R*(-(kC*rhoC) + kD*rhoD) + 
           151*bbCi^12*j2*pCi^4*R*(-(kC*rhoC) + kD*rhoD) + 
           bbCi^4*(4*kC*R*rhoC - 4*kD*R*rhoD) + 480*bbCi^7*dbbi*j2^2*pCi^2*
            xCi*yCi + 4*bbCi^6*(-(kC*pCi^2*R*rhoC) + 2*dpi*kD*pCi*R*rhoD + 
             kD*pCi^2*R*rhoD - 6*dyi*j2^2*xCi - 6*dxi*j2^2*yCi) + 
           8*bbCi^5*dbbi*(kD*pCi^2*R*rhoD - 6*j2^2*xCi*yCi) + 
           16*j2*R*(-(kC*rhoC) + kD*rhoD)*(xCi^2 - yCi^2) + 
           5*bbCi^8*j2*(-7*kC*R*rhoC + 7*kD*R*rhoD + 24*j2*pCi*
              (dyi*pCi*xCi + dxi*pCi*yCi + 2*dpi*xCi*yCi)))*Cos[2*tt] - 
         12*bbCi^12*dyi*j2^3*Cos[3*tt] + 144*bbCi^14*dyi*j2^3*pCi^2*
          Cos[3*tt] - 420*bbCi^16*dyi*j2^3*pCi^4*Cos[3*tt] + 
         112*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[3*tt] - 16*bbCi^8*j2^2*kC*pCi^2*R*
          rhoC*xCi*Cos[3*tt] - 112*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[3*tt] + 
         16*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[3*tt] - 72*bbCi^11*dbbi*j2^3*
          yCi*Cos[3*tt] + 288*bbCi^14*dpi*j2^3*pCi*yCi*Cos[3*tt] + 
         1152*bbCi^13*dbbi*j2^3*pCi^2*yCi*Cos[3*tt] - 1680*bbCi^16*dpi*j2^3*
          pCi^3*yCi*Cos[3*tt] - 4200*bbCi^15*dbbi*j2^3*pCi^4*yCi*Cos[3*tt] + 
         840*bbCi^10*j2^2*kC*R*rhoC*Cos[4*tt] - 1872*bbCi^12*j2^2*kC*pCi^2*R*
          rhoC*Cos[4*tt] + 1032*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[4*tt] - 
         840*bbCi^10*j2^2*kD*R*rhoD*Cos[4*tt] + 1872*bbCi^12*j2^2*kD*pCi^2*R*
          rhoD*Cos[4*tt] - 1032*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[4*tt] + 
         498*bbCi^10*j2^2*kC*R*rhoC*Cos[tt - 5*tti] - 1080*bbCi^12*j2^2*kC*
          pCi^2*R*rhoC*Cos[tt - 5*tti] + 582*bbCi^14*j2^2*kC*pCi^4*R*rhoC*
          Cos[tt - 5*tti] - 498*bbCi^10*j2^2*kD*R*rhoD*Cos[tt - 5*tti] + 
         1080*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[tt - 5*tti] - 
         582*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[tt - 5*tti] + 
         108*bbCi^12*dyi*j2^3*Cos[tt - 4*tti] - 648*bbCi^14*dyi*j2^3*pCi^2*
          Cos[tt - 4*tti] + 540*bbCi^16*dyi*j2^3*pCi^4*Cos[tt - 4*tti] - 
         144*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[tt - 4*tti] + 
         144*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[tt - 4*tti] + 
         144*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[tt - 4*tti] - 
         144*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[tt - 4*tti] + 
         648*bbCi^11*dbbi*j2^3*yCi*Cos[tt - 4*tti] - 1296*bbCi^14*dpi*j2^3*
          pCi*yCi*Cos[tt - 4*tti] - 5184*bbCi^13*dbbi*j2^3*pCi^2*yCi*
          Cos[tt - 4*tti] + 2160*bbCi^16*dpi*j2^3*pCi^3*yCi*Cos[tt - 4*tti] + 
         5400*bbCi^15*dbbi*j2^3*pCi^4*yCi*Cos[tt - 4*tti] - 
         336*bbCi^6*j2*kC*R*rhoC*Cos[tt - 3*tti] - 244*bbCi^10*j2^2*kC*R*rhoC*
          Cos[tt - 3*tti] + 336*bbCi^8*j2*kC*pCi^2*R*rhoC*Cos[tt - 3*tti] + 
         1608*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[tt - 3*tti] - 
         1508*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[tt - 3*tti] + 
         336*bbCi^6*j2*kD*R*rhoD*Cos[tt - 3*tti] + 244*bbCi^10*j2^2*kD*R*rhoD*
          Cos[tt - 3*tti] - 672*bbCi^8*dpi*j2*kD*pCi*R*rhoD*Cos[tt - 3*tti] - 
         336*bbCi^8*j2*kD*pCi^2*R*rhoD*Cos[tt - 3*tti] - 
         672*bbCi^7*dbbi*j2*kD*pCi^2*R*rhoD*Cos[tt - 3*tti] - 
         1608*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[tt - 3*tti] + 
         1508*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[tt - 3*tti] - 
         196*bbCi^10*j2^2*kC*R*rhoC*Cos[2*(tt - 3*tti)] + 
         392*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[2*(tt - 3*tti)] - 
         196*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[2*(tt - 3*tti)] + 
         196*bbCi^10*j2^2*kD*R*rhoD*Cos[2*(tt - 3*tti)] - 
         392*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[2*(tt - 3*tti)] + 
         196*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[2*(tt - 3*tti)] + 
         336*bbCi^12*dyi*j2^3*Cos[2*tt - 3*tti] - 2016*bbCi^14*dyi*j2^3*pCi^2*
          Cos[2*tt - 3*tti] + 1680*bbCi^16*dyi*j2^3*pCi^4*Cos[2*tt - 3*tti] + 
         448*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[2*tt - 3*tti] - 
         448*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[2*tt - 3*tti] - 
         448*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[2*tt - 3*tti] + 
         448*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[2*tt - 3*tti] + 
         2016*bbCi^11*dbbi*j2^3*yCi*Cos[2*tt - 3*tti] - 4032*bbCi^14*dpi*j2^3*
          pCi*yCi*Cos[2*tt - 3*tti] - 16128*bbCi^13*dbbi*j2^3*pCi^2*yCi*
          Cos[2*tt - 3*tti] + 6720*bbCi^16*dpi*j2^3*pCi^3*yCi*
          Cos[2*tt - 3*tti] + 16800*bbCi^15*dbbi*j2^3*pCi^4*yCi*
          Cos[2*tt - 3*tti] - 216*bbCi^12*dyi*j2^3*Cos[tt - 2*tti] + 
         792*bbCi^14*dyi*j2^3*pCi^2*Cos[tt - 2*tti] - 1440*bbCi^16*dyi*j2^3*
          pCi^4*Cos[tt - 2*tti] - 768*bbCi^6*j2^2*kC*R*rhoC*xCi*
          Cos[tt - 2*tti] + 672*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*
          Cos[tt - 2*tti] + 768*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[tt - 2*tti] - 
         672*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[tt - 2*tti] - 
         1296*bbCi^11*dbbi*j2^3*yCi*Cos[tt - 2*tti] + 1584*bbCi^14*dpi*j2^3*
          pCi*yCi*Cos[tt - 2*tti] + 6336*bbCi^13*dbbi*j2^3*pCi^2*yCi*
          Cos[tt - 2*tti] - 5760*bbCi^16*dpi*j2^3*pCi^3*yCi*Cos[tt - 2*tti] - 
         14400*bbCi^15*dbbi*j2^3*pCi^4*yCi*Cos[tt - 2*tti] + 
         672*bbCi^10*j2^2*kC*R*rhoC*Cos[2*(tt - 2*tti)] - 
         3360*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[2*(tt - 2*tti)] + 
         2688*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[2*(tt - 2*tti)] - 
         672*bbCi^10*j2^2*kD*R*rhoD*Cos[2*(tt - 2*tti)] + 
         3360*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[2*(tt - 2*tti)] - 
         2688*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[2*(tt - 2*tti)] + 
         480*bbCi^6*j2*kC*R*rhoC*Cos[tt - tti] + 762*bbCi^10*j2^2*kC*R*rhoC*
          Cos[tt - tti] - 1824*bbCi^8*j2*kC*pCi^2*R*rhoC*Cos[tt - tti] - 
         2688*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[tt - tti] + 
         3654*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[tt - tti] - 
         480*bbCi^6*j2*kD*R*rhoD*Cos[tt - tti] - 762*bbCi^10*j2^2*kD*R*rhoD*
          Cos[tt - tti] + 3648*bbCi^8*dpi*j2*kD*pCi*R*rhoD*Cos[tt - tti] + 
         1824*bbCi^8*j2*kD*pCi^2*R*rhoD*Cos[tt - tti] + 3648*bbCi^7*dbbi*j2*
          kD*pCi^2*R*rhoD*Cos[tt - tti] + 2688*bbCi^12*j2^2*kD*pCi^2*R*rhoD*
          Cos[tt - tti] - 3654*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[tt - tti] - 
         1200*bbCi^10*j2^2*kC*R*rhoC*Cos[2*(tt - tti)] + 
         4896*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[2*(tt - tti)] - 
         6576*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[2*(tt - tti)] + 
         1200*bbCi^10*j2^2*kD*R*rhoD*Cos[2*(tt - tti)] - 
         4896*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[2*(tt - tti)] + 
         6576*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[2*(tt - tti)] - 
         98*bbCi^10*j2^2*kC*R*rhoC*Cos[3*(tt - tti)] + 112*bbCi^12*j2^2*kC*
          pCi^2*R*rhoC*Cos[3*(tt - tti)] - 14*bbCi^14*j2^2*kC*pCi^4*R*rhoC*
          Cos[3*(tt - tti)] + 98*bbCi^10*j2^2*kD*R*rhoD*Cos[3*(tt - tti)] - 
         112*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[3*(tt - tti)] + 
         14*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[3*(tt - tti)] - 
         288*bbCi^12*dyi*j2^3*Cos[2*tt - tti] + 2304*bbCi^14*dyi*j2^3*pCi^2*
          Cos[2*tt - tti] - 4320*bbCi^16*dyi*j2^3*pCi^4*Cos[2*tt - tti] - 
         768*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[2*tt - tti] + 
         3072*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[2*tt - tti] + 
         768*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[2*tt - tti] - 
         3072*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[2*tt - tti] - 
         1728*bbCi^11*dbbi*j2^3*yCi*Cos[2*tt - tti] + 4608*bbCi^14*dpi*j2^3*
          pCi*yCi*Cos[2*tt - tti] + 18432*bbCi^13*dbbi*j2^3*pCi^2*yCi*
          Cos[2*tt - tti] - 17280*bbCi^16*dpi*j2^3*pCi^3*yCi*
          Cos[2*tt - tti] - 43200*bbCi^15*dbbi*j2^3*pCi^4*yCi*
          Cos[2*tt - tti] - 12*bbCi^10*j2^2*kC*R*rhoC*Cos[3*tt - tti] + 
         312*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[3*tt - tti] - 
         444*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[3*tt - tti] + 
         12*bbCi^10*j2^2*kD*R*rhoD*Cos[3*tt - tti] - 312*bbCi^12*j2^2*kD*
          pCi^2*R*rhoD*Cos[3*tt - tti] + 444*bbCi^14*j2^2*kD*pCi^4*R*rhoD*
          Cos[3*tt - tti] + 384*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[tti] + 
         384*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[tti] - 384*bbCi^6*j2^2*kD*R*
          rhoD*xCi*Cos[tti] - 384*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[tti] + 
         768*bbCi^6*j2*kC*R*rhoC*Cos[2*tti] + 480*bbCi^10*j2^2*kC*R*rhoC*
          Cos[2*tti] - 768*bbCi^8*j2*kC*pCi^2*R*rhoC*Cos[2*tti] - 
         2784*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[2*tti] + 2304*bbCi^14*j2^2*kC*
          pCi^4*R*rhoC*Cos[2*tti] - 768*bbCi^6*j2*kD*R*rhoD*Cos[2*tti] - 
         480*bbCi^10*j2^2*kD*R*rhoD*Cos[2*tti] + 1536*bbCi^8*dpi*j2*kD*pCi*R*
          rhoD*Cos[2*tti] + 768*bbCi^8*j2*kD*pCi^2*R*rhoD*Cos[2*tti] + 
         1536*bbCi^7*dbbi*j2*kD*pCi^2*R*rhoD*Cos[2*tti] + 
         2784*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[2*tti] - 2304*bbCi^14*j2^2*kD*
          pCi^4*R*rhoD*Cos[2*tti] + 640*bbCi^6*j2^2*kC*R*rhoC*xCi*
          Cos[3*tti] - 640*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[3*tti] - 
         640*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[3*tti] + 640*bbCi^8*j2^2*kD*pCi^2*
          R*rhoD*xCi*Cos[3*tti] - 672*bbCi^10*j2^2*kC*R*rhoC*Cos[4*tti] + 
         1344*bbCi^12*j2^2*kC*pCi^2*R*rhoC*Cos[4*tti] - 672*bbCi^14*j2^2*kC*
          pCi^4*R*rhoC*Cos[4*tti] + 672*bbCi^10*j2^2*kD*R*rhoD*Cos[4*tti] - 
         1344*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[4*tti] + 672*bbCi^14*j2^2*kD*
          pCi^4*R*rhoD*Cos[4*tti] - 144*bbCi^6*j2*kC*R*rhoC*Cos[tt + tti] - 
         336*bbCi^10*j2^2*kC*R*rhoC*Cos[tt + tti] + 144*bbCi^8*j2*kC*pCi^2*R*
          rhoC*Cos[tt + tti] + 1248*bbCi^12*j2^2*kC*pCi^2*R*rhoC*
          Cos[tt + tti] - 624*bbCi^14*j2^2*kC*pCi^4*R*rhoC*Cos[tt + tti] + 
         144*bbCi^6*j2*kD*R*rhoD*Cos[tt + tti] + 336*bbCi^10*j2^2*kD*R*rhoD*
          Cos[tt + tti] - 288*bbCi^8*dpi*j2*kD*pCi*R*rhoD*Cos[tt + tti] - 
         144*bbCi^8*j2*kD*pCi^2*R*rhoD*Cos[tt + tti] - 288*bbCi^7*dbbi*j2*kD*
          pCi^2*R*rhoD*Cos[tt + tti] - 1248*bbCi^12*j2^2*kD*pCi^2*R*rhoD*
          Cos[tt + tti] + 624*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[tt + tti] - 
         636*bbCi^10*j2^2*kC*R*rhoC*Cos[2*(tt + tti)] + 1464*bbCi^12*j2^2*kC*
          pCi^2*R*rhoC*Cos[2*(tt + tti)] - 828*bbCi^14*j2^2*kC*pCi^4*R*rhoC*
          Cos[2*(tt + tti)] + 636*bbCi^10*j2^2*kD*R*rhoD*Cos[2*(tt + tti)] - 
         1464*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[2*(tt + tti)] + 
         828*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[2*(tt + tti)] + 
         144*bbCi^12*dyi*j2^3*Cos[2*tt + tti] - 864*bbCi^14*dyi*j2^3*pCi^2*
          Cos[2*tt + tti] + 720*bbCi^16*dyi*j2^3*pCi^4*Cos[2*tt + tti] + 
         192*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[2*tt + tti] - 
         192*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[2*tt + tti] - 
         192*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[2*tt + tti] + 
         192*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[2*tt + tti] + 
         864*bbCi^11*dbbi*j2^3*yCi*Cos[2*tt + tti] - 1728*bbCi^14*dpi*j2^3*
          pCi*yCi*Cos[2*tt + tti] - 6912*bbCi^13*dbbi*j2^3*pCi^2*yCi*
          Cos[2*tt + tti] + 2880*bbCi^16*dpi*j2^3*pCi^3*yCi*Cos[2*tt + tti] + 
         7200*bbCi^15*dbbi*j2^3*pCi^4*yCi*Cos[2*tt + tti] - 
         42*bbCi^10*j2^2*kC*R*rhoC*Cos[3*tt + tti] + 48*bbCi^12*j2^2*kC*pCi^2*
          R*rhoC*Cos[3*tt + tti] - 6*bbCi^14*j2^2*kC*pCi^4*R*rhoC*
          Cos[3*tt + tti] + 42*bbCi^10*j2^2*kD*R*rhoD*Cos[3*tt + tti] - 
         48*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[3*tt + tti] + 
         6*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[3*tt + tti] + 
         72*bbCi^12*dyi*j2^3*Cos[tt + 2*tti] - 1224*bbCi^14*dyi*j2^3*pCi^2*
          Cos[tt + 2*tti] + 1440*bbCi^16*dyi*j2^3*pCi^4*Cos[tt + 2*tti] - 
         576*bbCi^6*j2^2*kC*R*rhoC*xCi*Cos[tt + 2*tti] + 
         480*bbCi^8*j2^2*kC*pCi^2*R*rhoC*xCi*Cos[tt + 2*tti] + 
         576*bbCi^6*j2^2*kD*R*rhoD*xCi*Cos[tt + 2*tti] - 
         480*bbCi^8*j2^2*kD*pCi^2*R*rhoD*xCi*Cos[tt + 2*tti] + 
         432*bbCi^11*dbbi*j2^3*yCi*Cos[tt + 2*tti] - 2448*bbCi^14*dpi*j2^3*
          pCi*yCi*Cos[tt + 2*tti] - 9792*bbCi^13*dbbi*j2^3*pCi^2*yCi*
          Cos[tt + 2*tti] + 5760*bbCi^16*dpi*j2^3*pCi^3*yCi*Cos[tt + 2*tti] + 
         14400*bbCi^15*dbbi*j2^3*pCi^4*yCi*Cos[tt + 2*tti] + 
         208*bbCi^10*j2^2*kC*R*rhoC*Cos[tt + 3*tti] - 296*bbCi^12*j2^2*kC*
          pCi^2*R*rhoC*Cos[tt + 3*tti] + 88*bbCi^14*j2^2*kC*pCi^4*R*rhoC*
          Cos[tt + 3*tti] - 208*bbCi^10*j2^2*kD*R*rhoD*Cos[tt + 3*tti] + 
         296*bbCi^12*j2^2*kD*pCi^2*R*rhoD*Cos[tt + 3*tti] - 
         88*bbCi^14*j2^2*kD*pCi^4*R*rhoD*Cos[tt + 3*tti] - 
         192*bbCi^8*dxi*j2^2*Sin[tt] - 384*bbCi^7*dbbi*dxi*j2^2*Sin[tt] - 
         276*bbCi^12*dxi*j2^3*Sin[tt] + 1920*bbCi^10*dpi*dxi*j2^2*pCi*
          Sin[tt] + 960*bbCi^10*dxi*j2^2*pCi^2*Sin[tt] + 
         3840*bbCi^9*dbbi*dxi*j2^2*pCi^2*Sin[tt] + 1776*bbCi^14*dxi*j2^3*
          pCi^2*Sin[tt] - 2940*bbCi^16*dxi*j2^3*pCi^4*Sin[tt] - 
         384*bbCi^2*dyi*j2*kD*R*rhoD*Sin[tt] - 384*bbCi^7*dbbi*j2^2*xCi*
          Sin[tt] - 192*bbCi^6*dbbi^2*j2^2*xCi*Sin[tt] + 
         960*bbCi^10*dpi^2*j2^2*xCi*Sin[tt] - 1656*bbCi^11*dbbi*j2^3*xCi*
          Sin[tt] + 1920*bbCi^10*dpi*j2^2*pCi*xCi*Sin[tt] + 
         7680*bbCi^9*dbbi*dpi*j2^2*pCi*xCi*Sin[tt] + 3552*bbCi^14*dpi*j2^3*
          pCi*xCi*Sin[tt] + 3840*bbCi^9*dbbi*j2^2*pCi^2*xCi*Sin[tt] + 
         5760*bbCi^8*dbbi^2*j2^2*pCi^2*xCi*Sin[tt] + 14208*bbCi^13*dbbi*j2^3*
          pCi^2*xCi*Sin[tt] - 11760*bbCi^16*dpi*j2^3*pCi^3*xCi*Sin[tt] - 
         29400*bbCi^15*dbbi*j2^3*pCi^4*xCi*Sin[tt] + 384*bbCi^2*j2*kC*R*rhoC*
          yCi*Sin[tt] + 32*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[tt] - 
         128*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[tt] - 384*bbCi^2*j2*kD*R*
          rhoD*yCi*Sin[tt] + 1536*bbCi*dbbi*j2*kD*R*rhoD*yCi*Sin[tt] - 
         32*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[tt] + 128*bbCi^8*j2^2*kD*pCi^2*R*
          rhoD*yCi*Sin[tt] - 1080*bbCi^15*dbbi*j2^3*Sin[2*tt] + 
         1944*bbCi^18*dpi*j2^3*pCi*Sin[2*tt] + 11664*bbCi^17*dbbi*j2^3*pCi^2*
          Sin[2*tt] - 9936*bbCi^20*dpi*j2^3*pCi^3*Sin[2*tt] - 
         34776*bbCi^19*dbbi*j2^3*pCi^4*Sin[2*tt] + 9720*bbCi^22*dpi*j2^3*
          pCi^5*Sin[2*tt] + 25920*bbCi^21*dbbi*j2^3*pCi^6*Sin[2*tt] + 
         384*bbCi^8*dxi*j2^3*xCi*Sin[2*tt] - 1920*bbCi^10*dxi*j2^3*pCi^2*xCi*
          Sin[2*tt] + 384*bbCi^7*dbbi*j2^3*xCi^2*Sin[2*tt] - 
         1920*bbCi^10*dpi*j2^3*pCi*xCi^2*Sin[2*tt] - 3840*bbCi^9*dbbi*j2^3*
          pCi^2*xCi^2*Sin[2*tt] - 384*bbCi^8*dyi*j2^3*yCi*Sin[2*tt] + 
         1920*bbCi^10*dyi*j2^3*pCi^2*yCi*Sin[2*tt] - 512*bbCi^2*j2^2*kC*R*
          rhoC*xCi*yCi*Sin[2*tt] + 512*bbCi^2*j2^2*kD*R*rhoD*xCi*yCi*
          Sin[2*tt] - 384*bbCi^7*dbbi*j2^3*yCi^2*Sin[2*tt] + 
         1920*bbCi^10*dpi*j2^3*pCi*yCi^2*Sin[2*tt] + 3840*bbCi^9*dbbi*j2^3*
          pCi^2*yCi^2*Sin[2*tt] + 12*bbCi^12*dxi*j2^3*Sin[3*tt] - 
         144*bbCi^14*dxi*j2^3*pCi^2*Sin[3*tt] + 420*bbCi^16*dxi*j2^3*pCi^4*
          Sin[3*tt] + 72*bbCi^11*dbbi*j2^3*xCi*Sin[3*tt] - 
         288*bbCi^14*dpi*j2^3*pCi*xCi*Sin[3*tt] - 1152*bbCi^13*dbbi*j2^3*
          pCi^2*xCi*Sin[3*tt] + 1680*bbCi^16*dpi*j2^3*pCi^3*xCi*Sin[3*tt] + 
         4200*bbCi^15*dbbi*j2^3*pCi^4*xCi*Sin[3*tt] + 112*bbCi^6*j2^2*kC*R*
          rhoC*yCi*Sin[3*tt] - 16*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[3*tt] - 
         112*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[3*tt] + 16*bbCi^8*j2^2*kD*pCi^2*R*
          rhoD*yCi*Sin[3*tt] + 675*bbCi^15*dbbi*j2^3*Sin[tt - 5*tti] - 
         2079*bbCi^18*dpi*j2^3*pCi*Sin[tt - 5*tti] - 12474*bbCi^17*dbbi*j2^3*
          pCi^2*Sin[tt - 5*tti] + 6498*bbCi^20*dpi*j2^3*pCi^3*
          Sin[tt - 5*tti] + 22743*bbCi^19*dbbi*j2^3*pCi^4*Sin[tt - 5*tti] - 
         3915*bbCi^22*dpi*j2^3*pCi^5*Sin[tt - 5*tti] - 10440*bbCi^21*dbbi*
          j2^3*pCi^6*Sin[tt - 5*tti] + 108*bbCi^12*dxi*j2^3*Sin[tt - 4*tti] - 
         648*bbCi^14*dxi*j2^3*pCi^2*Sin[tt - 4*tti] + 540*bbCi^16*dxi*j2^3*
          pCi^4*Sin[tt - 4*tti] + 648*bbCi^11*dbbi*j2^3*xCi*Sin[tt - 4*tti] - 
         1296*bbCi^14*dpi*j2^3*pCi*xCi*Sin[tt - 4*tti] - 
         5184*bbCi^13*dbbi*j2^3*pCi^2*xCi*Sin[tt - 4*tti] + 
         2160*bbCi^16*dpi*j2^3*pCi^3*xCi*Sin[tt - 4*tti] + 
         5400*bbCi^15*dbbi*j2^3*pCi^4*xCi*Sin[tt - 4*tti] + 
         144*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[tt - 4*tti] - 
         144*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[tt - 4*tti] - 
         144*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[tt - 4*tti] + 
         144*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[tt - 4*tti] + 
         1008*bbCi^11*dbbi*j2^2*Sin[tt - 3*tti] + 2520*bbCi^10*dbbi^2*j2^2*
          Sin[tt - 3*tti] - 1008*bbCi^14*dpi^2*j2^2*Sin[tt - 3*tti] + 
         1275*bbCi^15*dbbi*j2^3*Sin[tt - 3*tti] - 2016*bbCi^14*dpi*j2^2*pCi*
          Sin[tt - 3*tti] - 16128*bbCi^13*dbbi*dpi*j2^2*pCi*Sin[tt - 3*tti] - 
         513*bbCi^18*dpi*j2^3*pCi*Sin[tt - 3*tti] - 8064*bbCi^13*dbbi*j2^2*
          pCi^2*Sin[tt - 3*tti] - 28224*bbCi^12*dbbi^2*j2^2*pCi^2*
          Sin[tt - 3*tti] + 5040*bbCi^16*dpi^2*j2^2*pCi^2*Sin[tt - 3*tti] - 
         3078*bbCi^17*dbbi*j2^3*pCi^2*Sin[tt - 3*tti] + 3360*bbCi^16*dpi*j2^2*
          pCi^3*Sin[tt - 3*tti] + 33600*bbCi^15*dbbi*dpi*j2^2*pCi^3*
          Sin[tt - 3*tti] - 4902*bbCi^20*dpi*j2^3*pCi^3*Sin[tt - 3*tti] + 
         8400*bbCi^15*dbbi*j2^2*pCi^4*Sin[tt - 3*tti] + 37800*bbCi^14*dbbi^2*
          j2^2*pCi^4*Sin[tt - 3*tti] - 17157*bbCi^19*dbbi*j2^3*pCi^4*
          Sin[tt - 3*tti] + 5535*bbCi^22*dpi*j2^3*pCi^5*Sin[tt - 3*tti] + 
         14760*bbCi^21*dbbi*j2^3*pCi^6*Sin[tt - 3*tti] + 
         1470*bbCi^15*dbbi*j2^3*Sin[2*(tt - 3*tti)] - 2058*bbCi^18*dpi*j2^3*
          pCi*Sin[2*(tt - 3*tti)] - 12348*bbCi^17*dbbi*j2^3*pCi^2*
          Sin[2*(tt - 3*tti)] + 6468*bbCi^20*dpi*j2^3*pCi^3*
          Sin[2*(tt - 3*tti)] + 22638*bbCi^19*dbbi*j2^3*pCi^4*
          Sin[2*(tt - 3*tti)] - 4410*bbCi^22*dpi*j2^3*pCi^5*
          Sin[2*(tt - 3*tti)] - 11760*bbCi^21*dbbi*j2^3*pCi^6*
          Sin[2*(tt - 3*tti)] - 336*bbCi^12*dxi*j2^3*Sin[2*tt - 3*tti] + 
         2016*bbCi^14*dxi*j2^3*pCi^2*Sin[2*tt - 3*tti] - 
         1680*bbCi^16*dxi*j2^3*pCi^4*Sin[2*tt - 3*tti] - 
         2016*bbCi^11*dbbi*j2^3*xCi*Sin[2*tt - 3*tti] + 4032*bbCi^14*dpi*j2^3*
          pCi*xCi*Sin[2*tt - 3*tti] + 16128*bbCi^13*dbbi*j2^3*pCi^2*xCi*
          Sin[2*tt - 3*tti] - 6720*bbCi^16*dpi*j2^3*pCi^3*xCi*
          Sin[2*tt - 3*tti] - 16800*bbCi^15*dbbi*j2^3*pCi^4*xCi*
          Sin[2*tt - 3*tti] + 448*bbCi^6*j2^2*kC*R*rhoC*yCi*
          Sin[2*tt - 3*tti] - 448*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*
          Sin[2*tt - 3*tti] - 448*bbCi^6*j2^2*kD*R*rhoD*yCi*
          Sin[2*tt - 3*tti] + 448*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*
          Sin[2*tt - 3*tti] + 72*bbCi^12*dxi*j2^3*Sin[tt - 2*tti] + 
         360*bbCi^14*dxi*j2^3*pCi^2*Sin[tt - 2*tti] - 720*bbCi^16*dxi*j2^3*
          pCi^4*Sin[tt - 2*tti] + 432*bbCi^11*dbbi*j2^3*xCi*Sin[tt - 2*tti] + 
         720*bbCi^14*dpi*j2^3*pCi*xCi*Sin[tt - 2*tti] + 
         2880*bbCi^13*dbbi*j2^3*pCi^2*xCi*Sin[tt - 2*tti] - 
         2880*bbCi^16*dpi*j2^3*pCi^3*xCi*Sin[tt - 2*tti] - 
         7200*bbCi^15*dbbi*j2^3*pCi^4*xCi*Sin[tt - 2*tti] - 
         960*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[tt - 2*tti] + 
         1248*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[tt - 2*tti] + 
         960*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[tt - 2*tti] - 
         1248*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[tt - 2*tti] - 
         2520*bbCi^15*dbbi*j2^3*Sin[2*(tt - 2*tti)] + 4536*bbCi^18*dpi*j2^3*
          pCi*Sin[2*(tt - 2*tti)] + 27216*bbCi^17*dbbi*j2^3*pCi^2*
          Sin[2*(tt - 2*tti)] - 23184*bbCi^20*dpi*j2^3*pCi^3*
          Sin[2*(tt - 2*tti)] - 81144*bbCi^19*dbbi*j2^3*pCi^4*
          Sin[2*(tt - 2*tti)] + 22680*bbCi^22*dpi*j2^3*pCi^5*
          Sin[2*(tt - 2*tti)] + 60480*bbCi^21*dbbi*j2^3*pCi^6*
          Sin[2*(tt - 2*tti)] - 864*bbCi^11*dbbi*j2^2*Sin[tt - tti] - 
         2160*bbCi^10*dbbi^2*j2^2*Sin[tt - tti] + 1152*bbCi^14*dpi^2*j2^2*
          Sin[tt - tti] - 1125*bbCi^15*dbbi*j2^3*Sin[tt - tti] + 
         2304*bbCi^14*dpi*j2^2*pCi*Sin[tt - tti] + 18432*bbCi^13*dbbi*dpi*
          j2^2*pCi*Sin[tt - tti] + 1125*bbCi^18*dpi*j2^3*pCi*Sin[tt - tti] + 
         9216*bbCi^13*dbbi*j2^2*pCi^2*Sin[tt - tti] + 32256*bbCi^12*dbbi^2*
          j2^2*pCi^2*Sin[tt - tti] - 12960*bbCi^16*dpi^2*j2^2*pCi^2*
          Sin[tt - tti] + 6750*bbCi^17*dbbi*j2^3*pCi^2*Sin[tt - tti] - 
         8640*bbCi^16*dpi*j2^2*pCi^3*Sin[tt - tti] - 86400*bbCi^15*dbbi*dpi*
          j2^2*pCi^3*Sin[tt - tti] - 12510*bbCi^20*dpi*j2^3*pCi^3*
          Sin[tt - tti] - 21600*bbCi^15*dbbi*j2^2*pCi^4*Sin[tt - tti] - 
         97200*bbCi^14*dbbi^2*j2^2*pCi^4*Sin[tt - tti] - 
         43785*bbCi^19*dbbi*j2^3*pCi^4*Sin[tt - tti] + 16065*bbCi^22*dpi*j2^3*
          pCi^5*Sin[tt - tti] + 42840*bbCi^21*dbbi*j2^3*pCi^6*Sin[tt - tti] - 
         256*kC^2*R^2*rhoC^2*Sin[tt - tti] + 256*kD^2*R^2*rhoD^2*
          Sin[tt - tti] + 2340*bbCi^15*dbbi*j2^3*Sin[2*(tt - tti)] - 
         4140*bbCi^18*dpi*j2^3*pCi*Sin[2*(tt - tti)] - 24840*bbCi^17*dbbi*
          j2^3*pCi^2*Sin[2*(tt - tti)] + 22392*bbCi^20*dpi*j2^3*pCi^3*
          Sin[2*(tt - tti)] + 78372*bbCi^19*dbbi*j2^3*pCi^4*
          Sin[2*(tt - tti)] - 32940*bbCi^22*dpi*j2^3*pCi^5*
          Sin[2*(tt - tti)] - 87840*bbCi^21*dbbi*j2^3*pCi^6*
          Sin[2*(tt - tti)] - 105*bbCi^15*dbbi*j2^3*Sin[3*(tt - tti)] + 
         273*bbCi^18*dpi*j2^3*pCi*Sin[3*(tt - tti)] + 1638*bbCi^17*dbbi*j2^3*
          pCi^2*Sin[3*(tt - tti)] - 1974*bbCi^20*dpi*j2^3*pCi^3*
          Sin[3*(tt - tti)] - 6909*bbCi^19*dbbi*j2^3*pCi^4*
          Sin[3*(tt - tti)] + 2205*bbCi^22*dpi*j2^3*pCi^5*Sin[3*(tt - tti)] + 
         5880*bbCi^21*dbbi*j2^3*pCi^6*Sin[3*(tt - tti)] + 
         288*bbCi^12*dxi*j2^3*Sin[2*tt - tti] - 2304*bbCi^14*dxi*j2^3*pCi^2*
          Sin[2*tt - tti] + 4320*bbCi^16*dxi*j2^3*pCi^4*Sin[2*tt - tti] + 
         1728*bbCi^11*dbbi*j2^3*xCi*Sin[2*tt - tti] - 4608*bbCi^14*dpi*j2^3*
          pCi*xCi*Sin[2*tt - tti] - 18432*bbCi^13*dbbi*j2^3*pCi^2*xCi*
          Sin[2*tt - tti] + 17280*bbCi^16*dpi*j2^3*pCi^3*xCi*
          Sin[2*tt - tti] + 43200*bbCi^15*dbbi*j2^3*pCi^4*xCi*
          Sin[2*tt - tti] - 768*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[2*tt - tti] + 
         3072*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[2*tt - tti] + 
         768*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[2*tt - tti] - 
         3072*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[2*tt - tti] + 
         90*bbCi^15*dbbi*j2^3*Sin[3*tt - tti] - 270*bbCi^18*dpi*j2^3*pCi*
          Sin[3*tt - tti] - 1620*bbCi^17*dbbi*j2^3*pCi^2*Sin[3*tt - tti] + 
         2556*bbCi^20*dpi*j2^3*pCi^3*Sin[3*tt - tti] + 8946*bbCi^19*dbbi*j2^3*
          pCi^4*Sin[3*tt - tti] - 5670*bbCi^22*dpi*j2^3*pCi^5*
          Sin[3*tt - tti] - 15120*bbCi^21*dbbi*j2^3*pCi^6*Sin[3*tt - tti] - 
         1152*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[tti] + 1920*bbCi^8*j2^2*kC*pCi^2*
          R*rhoC*yCi*Sin[tti] + 1152*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[tti] - 
         1920*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[tti] + 
         640*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[3*tti] - 640*bbCi^8*j2^2*kC*pCi^2*
          R*rhoC*yCi*Sin[3*tti] - 640*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[3*tti] + 
         640*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[3*tti] + 
         432*bbCi^11*dbbi*j2^2*Sin[tt + tti] + 1080*bbCi^10*dbbi^2*j2^2*
          Sin[tt + tti] - 432*bbCi^14*dpi^2*j2^2*Sin[tt + tti] - 
         855*bbCi^15*dbbi*j2^3*Sin[tt + tti] - 864*bbCi^14*dpi*j2^2*pCi*
          Sin[tt + tti] - 6912*bbCi^13*dbbi*dpi*j2^2*pCi*Sin[tt + tti] + 
         3429*bbCi^18*dpi*j2^3*pCi*Sin[tt + tti] - 3456*bbCi^13*dbbi*j2^2*
          pCi^2*Sin[tt + tti] - 12096*bbCi^12*dbbi^2*j2^2*pCi^2*
          Sin[tt + tti] + 2160*bbCi^16*dpi^2*j2^2*pCi^2*Sin[tt + tti] + 
         20574*bbCi^17*dbbi*j2^3*pCi^2*Sin[tt + tti] + 1440*bbCi^16*dpi*j2^2*
          pCi^3*Sin[tt + tti] + 14400*bbCi^15*dbbi*dpi*j2^2*pCi^3*
          Sin[tt + tti] - 27666*bbCi^20*dpi*j2^3*pCi^3*Sin[tt + tti] + 
         3600*bbCi^15*dbbi*j2^2*pCi^4*Sin[tt + tti] + 16200*bbCi^14*dbbi^2*
          j2^2*pCi^4*Sin[tt + tti] - 96831*bbCi^19*dbbi*j2^3*pCi^4*
          Sin[tt + tti] + 31725*bbCi^22*dpi*j2^3*pCi^5*Sin[tt + tti] + 
         84600*bbCi^21*dbbi*j2^3*pCi^6*Sin[tt + tti] + 270*bbCi^15*dbbi*j2^3*
          Sin[2*(tt + tti)] - 378*bbCi^18*dpi*j2^3*pCi*Sin[2*(tt + tti)] - 
         2268*bbCi^17*dbbi*j2^3*pCi^2*Sin[2*(tt + tti)] + 
         1188*bbCi^20*dpi*j2^3*pCi^3*Sin[2*(tt + tti)] + 
         4158*bbCi^19*dbbi*j2^3*pCi^4*Sin[2*(tt + tti)] - 
         810*bbCi^22*dpi*j2^3*pCi^5*Sin[2*(tt + tti)] - 
         2160*bbCi^21*dbbi*j2^3*pCi^6*Sin[2*(tt + tti)] - 
         144*bbCi^12*dxi*j2^3*Sin[2*tt + tti] + 864*bbCi^14*dxi*j2^3*pCi^2*
          Sin[2*tt + tti] - 720*bbCi^16*dxi*j2^3*pCi^4*Sin[2*tt + tti] - 
         864*bbCi^11*dbbi*j2^3*xCi*Sin[2*tt + tti] + 1728*bbCi^14*dpi*j2^3*
          pCi*xCi*Sin[2*tt + tti] + 6912*bbCi^13*dbbi*j2^3*pCi^2*xCi*
          Sin[2*tt + tti] - 2880*bbCi^16*dpi*j2^3*pCi^3*xCi*Sin[2*tt + tti] - 
         7200*bbCi^15*dbbi*j2^3*pCi^4*xCi*Sin[2*tt + tti] + 
         192*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[2*tt + tti] - 
         192*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[2*tt + tti] - 
         192*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[2*tt + tti] + 
         192*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[2*tt + tti] - 
         45*bbCi^15*dbbi*j2^3*Sin[3*tt + tti] + 117*bbCi^18*dpi*j2^3*pCi*
          Sin[3*tt + tti] + 702*bbCi^17*dbbi*j2^3*pCi^2*Sin[3*tt + tti] - 
         846*bbCi^20*dpi*j2^3*pCi^3*Sin[3*tt + tti] - 2961*bbCi^19*dbbi*j2^3*
          pCi^4*Sin[3*tt + tti] + 945*bbCi^22*dpi*j2^3*pCi^5*
          Sin[3*tt + tti] + 2520*bbCi^21*dbbi*j2^3*pCi^6*Sin[3*tt + tti] - 
         72*bbCi^12*dxi*j2^3*Sin[tt + 2*tti] + 1224*bbCi^14*dxi*j2^3*pCi^2*
          Sin[tt + 2*tti] - 1440*bbCi^16*dxi*j2^3*pCi^4*Sin[tt + 2*tti] - 
         432*bbCi^11*dbbi*j2^3*xCi*Sin[tt + 2*tti] + 2448*bbCi^14*dpi*j2^3*
          pCi*xCi*Sin[tt + 2*tti] + 9792*bbCi^13*dbbi*j2^3*pCi^2*xCi*
          Sin[tt + 2*tti] - 5760*bbCi^16*dpi*j2^3*pCi^3*xCi*Sin[tt + 2*tti] - 
         14400*bbCi^15*dbbi*j2^3*pCi^4*xCi*Sin[tt + 2*tti] - 
         576*bbCi^6*j2^2*kC*R*rhoC*yCi*Sin[tt + 2*tti] + 
         480*bbCi^8*j2^2*kC*pCi^2*R*rhoC*yCi*Sin[tt + 2*tti] + 
         576*bbCi^6*j2^2*kD*R*rhoD*yCi*Sin[tt + 2*tti] - 
         480*bbCi^8*j2^2*kD*pCi^2*R*rhoD*yCi*Sin[tt + 2*tti] - 
         90*bbCi^15*dbbi*j2^3*Sin[tt + 3*tti] + 78*bbCi^18*dpi*j2^3*pCi*
          Sin[tt + 3*tti] + 468*bbCi^17*dbbi*j2^3*pCi^2*Sin[tt + 3*tti] - 
         780*bbCi^20*dpi*j2^3*pCi^3*Sin[tt + 3*tti] - 2730*bbCi^19*dbbi*j2^3*
          pCi^4*Sin[tt + 3*tti] + 990*bbCi^22*dpi*j2^3*pCi^5*
          Sin[tt + 3*tti] + 2640*bbCi^21*dbbi*j2^3*pCi^6*Sin[tt + 3*tti]) + 
       384*(tt - tti)^2*(-384*bbCi*kC^2*R^2*rhoC^2 - 96*bbCi^5*j2*kC^2*R^2*
          rhoC^2 + 1440*bbCi^7*j2*kC^2*pCi^2*R^2*rhoC^2 + 
         384*bbCi*kD^2*R^2*rhoD^2 - 2304*dbbi*kD^2*R^2*rhoD^2 + 
         96*bbCi^5*j2*kD^2*R^2*rhoD^2 - 1440*bbCi^7*j2*kD^2*pCi^2*R^2*
          rhoD^2 - 8*bbCi*j2*(3024*bbCi^17*dbbi*j2^3*pCi^2*xCi - 
           16380*bbCi^19*dbbi*j2^3*pCi^4*xCi + 28800*bbCi^21*dbbi*j2^3*pCi^6*
            xCi + 80*R^2*(-(kC^2*rhoC^2) + kD^2*rhoD^2)*xCi - 
           108*bbCi^11*dbbi*j2^2*(dxi + xCi) + 252*bbCi^18*j2^3*pCi*
            (dxi*pCi + 2*dpi*xCi) - 1170*bbCi^20*j2^3*pCi^3*
            (dxi*pCi + 4*dpi*xCi) + 1800*bbCi^22*j2^3*pCi^5*
            (dxi*pCi + 6*dpi*xCi) + 1440*bbCi^13*dbbi*j2^2*pCi*
            (dxi*pCi + (2*dpi + pCi)*xCi) - 18*bbCi^16*j2^2*
            (dxi*(j2 + 25*pCi^3*(4*dpi + pCi)) + 50*dpi*pCi^2*(3*dpi + 2*pCi)*
              xCi) - 180*bbCi^15*dbbi*j2^2*(25*dxi*pCi^4 + 
             (j2 + 25*pCi^3*(4*dpi + pCi))*xCi) + 480*bbCi^7*dbbi*j2*kD*pCi^2*
            R*rhoD*yCi - 48*bbCi^6*j2*R*(dyi*kD*rhoD - kC*rhoC*yCi + 
             kD*rhoD*yCi) - 45*bbCi^10*j2^2*(6*dbbi^2*xCi + 
             R*(-(kC*rhoC) + kD*rhoD)*yCi) + 240*bbCi^8*j2*pCi*R*
            (dyi*kD*pCi*rhoD + (-(kC*pCi*rhoC) + 2*dpi*kD*rhoD + kD*pCi*rhoD)*
              yCi) + 15*bbCi^14*j2^2*(12*dxi*pCi^2 + 12*dpi^2*xCi + 
             24*dpi*pCi*(dxi + xCi) + pCi^4*(-1350*dbbi^2*xCi + 13*R*
                (kC*rhoC - kD*rhoD)*yCi)) + 6*bbCi^12*j2^2*
            (-3*dxi + 4*pCi^2*(210*dbbi^2*xCi + 13*R*(-(kC*rhoC) + kD*rhoD)*
                yCi)))*Cos[tt] + 4*bbCi^5*j2*(567*bbCi^15*dbbi*j2^3 - 
           1134*bbCi^18*dpi*j2^3*pCi - 9072*bbCi^17*dbbi*j2^3*pCi^2 + 
           11016*bbCi^20*dpi*j2^3*pCi^3 + 49572*bbCi^19*dbbi*j2^3*pCi^4 - 
           31590*bbCi^22*dpi*j2^3*pCi^5 - 105300*bbCi^21*dbbi*j2^3*pCi^6 + 
           24300*bbCi^24*dpi*j2^3*pCi^7 + 66825*bbCi^23*dbbi*j2^3*pCi^8 + 
           232*R^2*(kC^2*rhoC^2 - kD^2*rhoD^2) + 1440*bbCi^4*j2^2*pCi^2*R*
            (-(kC*rhoC) + kD*rhoD)*xCi*yCi - 144*bbCi^8*j2^3*
            (dxi*xCi - dyi*yCi) + 8*bbCi^2*R*(-(kC*rhoC) + kD*rhoD)*
            (29*kC*pCi^2*R*rhoC + 29*kD*pCi^2*R*rhoD - 36*j2^2*xCi*yCi) - 
           432*bbCi^7*dbbi*j2^3*(xCi^2 - yCi^2) + 5760*bbCi^9*dbbi*j2^3*pCi^2*
            (xCi^2 - yCi^2) - 18000*bbCi^11*dbbi*j2^3*pCi^4*(xCi^2 - yCi^2) - 
           3600*bbCi^12*j2^3*pCi^3*(dxi*pCi*xCi + 2*dpi*xCi^2 - dyi*pCi*yCi - 
             2*dpi*yCi^2) + 1440*bbCi^10*j2^3*pCi*(dxi*pCi*xCi + dpi*xCi^2 - 
             dyi*pCi*yCi - dpi*yCi^2))*Cos[2*tt] - 72*bbCi^17*dxi*j2^4*
          Cos[3*tt] + 792*bbCi^19*dxi*j2^4*pCi^2*Cos[3*tt] - 
         2520*bbCi^21*dxi*j2^4*pCi^4*Cos[3*tt] + 1800*bbCi^23*dxi*j2^4*pCi^6*
          Cos[3*tt] - 720*bbCi^16*dbbi*j2^4*xCi*Cos[3*tt] + 
         1584*bbCi^19*dpi*j2^4*pCi*xCi*Cos[3*tt] + 9504*bbCi^18*dbbi*j2^4*
          pCi^2*xCi*Cos[3*tt] - 10080*bbCi^21*dpi*j2^4*pCi^3*xCi*Cos[3*tt] - 
         35280*bbCi^20*dbbi*j2^4*pCi^4*xCi*Cos[3*tt] + 10800*bbCi^23*dpi*j2^4*
          pCi^5*xCi*Cos[3*tt] + 28800*bbCi^22*dbbi*j2^4*pCi^6*xCi*Cos[3*tt] - 
         264*bbCi^11*j2^3*kC*R*rhoC*yCi*Cos[3*tt] + 1440*bbCi^13*j2^3*kC*
          pCi^2*R*rhoC*yCi*Cos[3*tt] - 600*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*
          Cos[3*tt] + 264*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[3*tt] - 
         1440*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[3*tt] + 
         600*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[3*tt] + 
         1323*bbCi^20*dbbi*j2^4*Cos[tt - 5*tti] - 2268*bbCi^23*dpi*j2^4*pCi*
          Cos[tt - 5*tti] - 18144*bbCi^22*dbbi*j2^4*pCi^2*Cos[tt - 5*tti] + 
         17388*bbCi^25*dpi*j2^4*pCi^3*Cos[tt - 5*tti] + 78246*bbCi^24*dbbi*
          j2^4*pCi^4*Cos[tt - 5*tti] - 34020*bbCi^27*dpi*j2^4*pCi^5*
          Cos[tt - 5*tti] - 113400*bbCi^26*dbbi*j2^4*pCi^6*Cos[tt - 5*tti] + 
         18900*bbCi^29*dpi*j2^4*pCi^7*Cos[tt - 5*tti] + 51975*bbCi^28*dbbi*
          j2^4*pCi^8*Cos[tt - 5*tti] - 216*bbCi^11*j2^3*kC*R*rhoC*yCi*
          Cos[tt - 4*tti] + 1296*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*
          Cos[tt - 4*tti] - 1080*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*
          Cos[tt - 4*tti] + 216*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[tt - 4*tti] - 
         1296*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[tt - 4*tti] + 
         1080*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[tt - 4*tti] - 
         1260*bbCi^16*dbbi*j2^3*Cos[tt - 3*tti] - 5670*bbCi^15*dbbi^2*j2^3*
          Cos[tt - 3*tti] + 1386*bbCi^19*dpi^2*j2^3*Cos[tt - 3*tti] - 
         3780*bbCi^20*dbbi*j2^4*Cos[tt - 3*tti] + 2772*bbCi^19*dpi*j2^3*pCi*
          Cos[tt - 3*tti] + 33264*bbCi^18*dbbi*dpi*j2^3*pCi*Cos[tt - 3*tti] + 
         7560*bbCi^23*dpi*j2^4*pCi*Cos[tt - 3*tti] + 16632*bbCi^18*dbbi*j2^3*
          pCi^2*Cos[tt - 3*tti] + 91476*bbCi^17*dbbi^2*j2^3*pCi^2*
          Cos[tt - 3*tti] - 26460*bbCi^21*dpi^2*j2^3*pCi^2*Cos[tt - 3*tti] + 
         60480*bbCi^22*dbbi*j2^4*pCi^2*Cos[tt - 3*tti] - 
         17640*bbCi^21*dpi*j2^3*pCi^3*Cos[tt - 3*tti] - 246960*bbCi^20*dbbi*
          dpi*j2^3*pCi^3*Cos[tt - 3*tti] - 73440*bbCi^25*dpi*j2^4*pCi^3*
          Cos[tt - 3*tti] - 61740*bbCi^20*dbbi*j2^3*pCi^4*Cos[tt - 3*tti] - 
         401310*bbCi^19*dbbi^2*j2^3*pCi^4*Cos[tt - 3*tti] + 
         47250*bbCi^23*dpi^2*j2^3*pCi^4*Cos[tt - 3*tti] - 
         330480*bbCi^24*dbbi*j2^4*pCi^4*Cos[tt - 3*tti] + 
         18900*bbCi^23*dpi*j2^3*pCi^5*Cos[tt - 3*tti] + 302400*bbCi^22*dbbi*
          dpi*j2^3*pCi^5*Cos[tt - 3*tti] + 210600*bbCi^27*dpi*j2^4*pCi^5*
          Cos[tt - 3*tti] + 50400*bbCi^22*dbbi*j2^3*pCi^6*Cos[tt - 3*tti] + 
         378000*bbCi^21*dbbi^2*j2^3*pCi^6*Cos[tt - 3*tti] + 
         702000*bbCi^26*dbbi*j2^4*pCi^6*Cos[tt - 3*tti] - 
         162000*bbCi^29*dpi*j2^4*pCi^7*Cos[tt - 3*tti] - 
         445500*bbCi^28*dbbi*j2^4*pCi^8*Cos[tt - 3*tti] - 
         560*bbCi^5*j2*kC^2*R^2*rhoC^2*Cos[tt - 3*tti] + 
         560*bbCi^7*j2*kC^2*pCi^2*R^2*rhoC^2*Cos[tt - 3*tti] + 
         560*bbCi^5*j2*kD^2*R^2*rhoD^2*Cos[tt - 3*tti] - 
         560*bbCi^7*j2*kD^2*pCi^2*R^2*rhoD^2*Cos[tt - 3*tti] - 
         3087*bbCi^20*dbbi*j2^4*Cos[2*(tt - 3*tti)] + 5292*bbCi^23*dpi*j2^4*
          pCi*Cos[2*(tt - 3*tti)] + 42336*bbCi^22*dbbi*j2^4*pCi^2*
          Cos[2*(tt - 3*tti)] - 40572*bbCi^25*dpi*j2^4*pCi^3*
          Cos[2*(tt - 3*tti)] - 182574*bbCi^24*dbbi*j2^4*pCi^4*
          Cos[2*(tt - 3*tti)] + 79380*bbCi^27*dpi*j2^4*pCi^5*
          Cos[2*(tt - 3*tti)] + 264600*bbCi^26*dbbi*j2^4*pCi^6*
          Cos[2*(tt - 3*tti)] - 44100*bbCi^29*dpi*j2^4*pCi^7*
          Cos[2*(tt - 3*tti)] - 121275*bbCi^28*dbbi*j2^4*pCi^8*
          Cos[2*(tt - 3*tti)] + 504*bbCi^17*dxi*j2^4*Cos[2*tt - 3*tti] - 
         5544*bbCi^19*dxi*j2^4*pCi^2*Cos[2*tt - 3*tti] + 
         17640*bbCi^21*dxi*j2^4*pCi^4*Cos[2*tt - 3*tti] - 
         12600*bbCi^23*dxi*j2^4*pCi^6*Cos[2*tt - 3*tti] + 
         5040*bbCi^16*dbbi*j2^4*xCi*Cos[2*tt - 3*tti] - 
         11088*bbCi^19*dpi*j2^4*pCi*xCi*Cos[2*tt - 3*tti] - 
         66528*bbCi^18*dbbi*j2^4*pCi^2*xCi*Cos[2*tt - 3*tti] + 
         70560*bbCi^21*dpi*j2^4*pCi^3*xCi*Cos[2*tt - 3*tti] + 
         246960*bbCi^20*dbbi*j2^4*pCi^4*xCi*Cos[2*tt - 3*tti] - 
         75600*bbCi^23*dpi*j2^4*pCi^5*xCi*Cos[2*tt - 3*tti] - 
         201600*bbCi^22*dbbi*j2^4*pCi^6*xCi*Cos[2*tt - 3*tti] - 
         1008*bbCi^11*j2^3*kC*R*rhoC*yCi*Cos[2*tt - 3*tti] + 
         6048*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*Cos[2*tt - 3*tti] - 
         5040*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*Cos[2*tt - 3*tti] + 
         1008*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[2*tt - 3*tti] - 
         6048*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[2*tt - 3*tti] + 
         5040*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[2*tt - 3*tti] - 
         108*bbCi^17*dxi*j2^4*Cos[tt - 2*tti] + 1188*bbCi^19*dxi*j2^4*pCi^2*
          Cos[tt - 2*tti] - 3780*bbCi^21*dxi*j2^4*pCi^4*Cos[tt - 2*tti] + 
         2700*bbCi^23*dxi*j2^4*pCi^6*Cos[tt - 2*tti] - 1080*bbCi^16*dbbi*j2^4*
          xCi*Cos[tt - 2*tti] + 2376*bbCi^19*dpi*j2^4*pCi*xCi*
          Cos[tt - 2*tti] + 14256*bbCi^18*dbbi*j2^4*pCi^2*xCi*
          Cos[tt - 2*tti] - 15120*bbCi^21*dpi*j2^4*pCi^3*xCi*
          Cos[tt - 2*tti] - 52920*bbCi^20*dbbi*j2^4*pCi^4*xCi*
          Cos[tt - 2*tti] + 16200*bbCi^23*dpi*j2^4*pCi^5*xCi*
          Cos[tt - 2*tti] + 43200*bbCi^22*dbbi*j2^4*pCi^6*xCi*
          Cos[tt - 2*tti] + 720*bbCi^11*j2^3*kC*R*rhoC*yCi*Cos[tt - 2*tti] - 
         3312*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*Cos[tt - 2*tti] + 
         4320*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*Cos[tt - 2*tti] - 
         720*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[tt - 2*tti] + 
         3312*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[tt - 2*tti] - 
         4320*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[tt - 2*tti] + 
         5292*bbCi^20*dbbi*j2^4*Cos[2*(tt - 2*tti)] - 10584*bbCi^23*dpi*j2^4*
          pCi*Cos[2*(tt - 2*tti)] - 84672*bbCi^22*dbbi*j2^4*pCi^2*
          Cos[2*(tt - 2*tti)] + 102816*bbCi^25*dpi*j2^4*pCi^3*
          Cos[2*(tt - 2*tti)] + 462672*bbCi^24*dbbi*j2^4*pCi^4*
          Cos[2*(tt - 2*tti)] - 294840*bbCi^27*dpi*j2^4*pCi^5*
          Cos[2*(tt - 2*tti)] - 982800*bbCi^26*dbbi*j2^4*pCi^6*
          Cos[2*(tt - 2*tti)] + 226800*bbCi^29*dpi*j2^4*pCi^7*
          Cos[2*(tt - 2*tti)] + 623700*bbCi^28*dbbi*j2^4*pCi^8*
          Cos[2*(tt - 2*tti)] + 1080*bbCi^16*dbbi*j2^3*Cos[tt - tti] + 
         4860*bbCi^15*dbbi^2*j2^3*Cos[tt - tti] - 1404*bbCi^19*dpi^2*j2^3*
          Cos[tt - tti] + 4536*bbCi^20*dbbi*j2^4*Cos[tt - tti] - 
         2808*bbCi^19*dpi*j2^3*pCi*Cos[tt - tti] - 33696*bbCi^18*dbbi*dpi*
          j2^3*pCi*Cos[tt - tti] - 9072*bbCi^23*dpi*j2^4*pCi*Cos[tt - tti] - 
         16848*bbCi^18*dbbi*j2^3*pCi^2*Cos[tt - tti] - 92664*bbCi^17*dbbi^2*
          j2^3*pCi^2*Cos[tt - tti] + 35640*bbCi^21*dpi^2*j2^3*pCi^2*
          Cos[tt - tti] - 72576*bbCi^22*dbbi*j2^4*pCi^2*Cos[tt - tti] + 
         23760*bbCi^21*dpi*j2^3*pCi^3*Cos[tt - tti] + 332640*bbCi^20*dbbi*dpi*
          j2^3*pCi^3*Cos[tt - tti] + 90720*bbCi^25*dpi*j2^4*pCi^3*
          Cos[tt - tti] + 83160*bbCi^20*dbbi*j2^3*pCi^4*Cos[tt - tti] + 
         540540*bbCi^19*dbbi^2*j2^3*pCi^4*Cos[tt - tti] - 
         121500*bbCi^23*dpi^2*j2^3*pCi^4*Cos[tt - tti] + 
         408240*bbCi^24*dbbi*j2^4*pCi^4*Cos[tt - tti] - 
         48600*bbCi^23*dpi*j2^3*pCi^5*Cos[tt - tti] - 777600*bbCi^22*dbbi*dpi*
          j2^3*pCi^5*Cos[tt - tti] - 291600*bbCi^27*dpi*j2^4*pCi^5*
          Cos[tt - tti] - 129600*bbCi^22*dbbi*j2^3*pCi^6*Cos[tt - tti] - 
         972000*bbCi^21*dbbi^2*j2^3*pCi^6*Cos[tt - tti] - 
         972000*bbCi^26*dbbi*j2^4*pCi^6*Cos[tt - tti] + 324000*bbCi^29*dpi*
          j2^4*pCi^7*Cos[tt - tti] + 891000*bbCi^28*dbbi*j2^4*pCi^8*
          Cos[tt - tti] + 864*bbCi^5*j2*kC^2*R^2*rhoC^2*Cos[tt - tti] - 
         3360*bbCi^7*j2*kC^2*pCi^2*R^2*rhoC^2*Cos[tt - tti] - 
         864*bbCi^5*j2*kD^2*R^2*rhoD^2*Cos[tt - tti] + 3360*bbCi^7*j2*kD^2*
          pCi^2*R^2*rhoD^2*Cos[tt - tti] - 4914*bbCi^20*dbbi*j2^4*
          Cos[2*(tt - tti)] + 9720*bbCi^23*dpi*j2^4*pCi*Cos[2*(tt - tti)] + 
         77760*bbCi^22*dbbi*j2^4*pCi^2*Cos[2*(tt - tti)] - 
         95688*bbCi^25*dpi*j2^4*pCi^3*Cos[2*(tt - tti)] - 
         430596*bbCi^24*dbbi*j2^4*pCi^4*Cos[2*(tt - tti)] + 
         301320*bbCi^27*dpi*j2^4*pCi^5*Cos[2*(tt - tti)] + 
         1004400*bbCi^26*dbbi*j2^4*pCi^6*Cos[2*(tt - tti)] - 
         329400*bbCi^29*dpi*j2^4*pCi^7*Cos[2*(tt - tti)] - 
         905850*bbCi^28*dbbi*j2^4*pCi^8*Cos[2*(tt - tti)] + 
         882*bbCi^20*dbbi*j2^4*Cos[3*(tt - tti)] - 1512*bbCi^23*dpi*j2^4*pCi*
          Cos[3*(tt - tti)] - 12096*bbCi^22*dbbi*j2^4*pCi^2*
          Cos[3*(tt - tti)] + 11592*bbCi^25*dpi*j2^4*pCi^3*
          Cos[3*(tt - tti)] + 52164*bbCi^24*dbbi*j2^4*pCi^4*
          Cos[3*(tt - tti)] - 22680*bbCi^27*dpi*j2^4*pCi^5*
          Cos[3*(tt - tti)] - 75600*bbCi^26*dbbi*j2^4*pCi^6*
          Cos[3*(tt - tti)] + 12600*bbCi^29*dpi*j2^4*pCi^7*
          Cos[3*(tt - tti)] + 34650*bbCi^28*dbbi*j2^4*pCi^8*
          Cos[3*(tt - tti)] - 432*bbCi^17*dxi*j2^4*Cos[2*tt - tti] + 
         5616*bbCi^19*dxi*j2^4*pCi^2*Cos[2*tt - tti] - 23760*bbCi^21*dxi*j2^4*
          pCi^4*Cos[2*tt - tti] + 32400*bbCi^23*dxi*j2^4*pCi^6*
          Cos[2*tt - tti] - 4320*bbCi^16*dbbi*j2^4*xCi*Cos[2*tt - tti] + 
         11232*bbCi^19*dpi*j2^4*pCi*xCi*Cos[2*tt - tti] + 
         67392*bbCi^18*dbbi*j2^4*pCi^2*xCi*Cos[2*tt - tti] - 
         95040*bbCi^21*dpi*j2^4*pCi^3*xCi*Cos[2*tt - tti] - 
         332640*bbCi^20*dbbi*j2^4*pCi^4*xCi*Cos[2*tt - tti] + 
         194400*bbCi^23*dpi*j2^4*pCi^5*xCi*Cos[2*tt - tti] + 
         518400*bbCi^22*dbbi*j2^4*pCi^6*xCi*Cos[2*tt - tti] + 
         1296*bbCi^11*j2^3*kC*R*rhoC*yCi*Cos[2*tt - tti] - 
         11232*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*Cos[2*tt - tti] + 
         23760*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*Cos[2*tt - tti] - 
         1296*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[2*tt - tti] + 
         11232*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[2*tt - tti] - 
         23760*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[2*tt - tti] - 
         756*bbCi^20*dbbi*j2^4*Cos[3*tt - tti] + 1512*bbCi^23*dpi*j2^4*pCi*
          Cos[3*tt - tti] + 12096*bbCi^22*dbbi*j2^4*pCi^2*Cos[3*tt - tti] - 
         14688*bbCi^25*dpi*j2^4*pCi^3*Cos[3*tt - tti] - 66096*bbCi^24*dbbi*
          j2^4*pCi^4*Cos[3*tt - tti] + 42120*bbCi^27*dpi*j2^4*pCi^5*
          Cos[3*tt - tti] + 140400*bbCi^26*dbbi*j2^4*pCi^6*Cos[3*tt - tti] - 
         32400*bbCi^29*dpi*j2^4*pCi^7*Cos[3*tt - tti] - 89100*bbCi^28*dbbi*
          j2^4*pCi^8*Cos[3*tt - tti] + 144*bbCi^11*j2^3*kC*R*rhoC*yCi*
          Cos[tti] - 1440*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*Cos[tti] + 
         3600*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*Cos[tti] - 
         144*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[tti] + 1440*bbCi^13*j2^3*kD*pCi^2*
          R*rhoD*yCi*Cos[tti] - 3600*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*
          Cos[tti] + 1728*bbCi^5*j2*kC^2*R^2*rhoC^2*Cos[2*tti] - 
         1728*bbCi^7*j2*kC^2*pCi^2*R^2*rhoC^2*Cos[2*tti] - 
         1728*bbCi^5*j2*kD^2*R^2*rhoD^2*Cos[2*tti] + 1728*bbCi^7*j2*kD^2*
          pCi^2*R^2*rhoD^2*Cos[2*tti] - 540*bbCi^16*dbbi*j2^3*Cos[tt + tti] - 
         2430*bbCi^15*dbbi^2*j2^3*Cos[tt + tti] + 594*bbCi^19*dpi^2*j2^3*
          Cos[tt + tti] - 3024*bbCi^20*dbbi*j2^4*Cos[tt + tti] + 
         1188*bbCi^19*dpi*j2^3*pCi*Cos[tt + tti] + 14256*bbCi^18*dbbi*dpi*
          j2^3*pCi*Cos[tt + tti] + 6048*bbCi^23*dpi*j2^4*pCi*Cos[tt + tti] + 
         7128*bbCi^18*dbbi*j2^3*pCi^2*Cos[tt + tti] + 39204*bbCi^17*dbbi^2*
          j2^3*pCi^2*Cos[tt + tti] - 11340*bbCi^21*dpi^2*j2^3*pCi^2*
          Cos[tt + tti] + 48384*bbCi^22*dbbi*j2^4*pCi^2*Cos[tt + tti] - 
         7560*bbCi^21*dpi*j2^3*pCi^3*Cos[tt + tti] - 105840*bbCi^20*dbbi*dpi*
          j2^3*pCi^3*Cos[tt + tti] - 58752*bbCi^25*dpi*j2^4*pCi^3*
          Cos[tt + tti] - 26460*bbCi^20*dbbi*j2^3*pCi^4*Cos[tt + tti] - 
         171990*bbCi^19*dbbi^2*j2^3*pCi^4*Cos[tt + tti] + 
         20250*bbCi^23*dpi^2*j2^3*pCi^4*Cos[tt + tti] - 264384*bbCi^24*dbbi*
          j2^4*pCi^4*Cos[tt + tti] + 8100*bbCi^23*dpi*j2^3*pCi^5*
          Cos[tt + tti] + 129600*bbCi^22*dbbi*dpi*j2^3*pCi^5*Cos[tt + tti] + 
         168480*bbCi^27*dpi*j2^4*pCi^5*Cos[tt + tti] + 21600*bbCi^22*dbbi*
          j2^3*pCi^6*Cos[tt + tti] + 162000*bbCi^21*dbbi^2*j2^3*pCi^6*
          Cos[tt + tti] + 561600*bbCi^26*dbbi*j2^4*pCi^6*Cos[tt + tti] - 
         129600*bbCi^29*dpi*j2^4*pCi^7*Cos[tt + tti] - 356400*bbCi^28*dbbi*
          j2^4*pCi^8*Cos[tt + tti] - 240*bbCi^5*j2*kC^2*R^2*rhoC^2*
          Cos[tt + tti] + 240*bbCi^7*j2*kC^2*pCi^2*R^2*rhoC^2*Cos[tt + tti] + 
         240*bbCi^5*j2*kD^2*R^2*rhoD^2*Cos[tt + tti] - 240*bbCi^7*j2*kD^2*
          pCi^2*R^2*rhoD^2*Cos[tt + tti] - 567*bbCi^20*dbbi*j2^4*
          Cos[2*(tt + tti)] + 972*bbCi^23*dpi*j2^4*pCi*Cos[2*(tt + tti)] + 
         7776*bbCi^22*dbbi*j2^4*pCi^2*Cos[2*(tt + tti)] - 
         7452*bbCi^25*dpi*j2^4*pCi^3*Cos[2*(tt + tti)] - 
         33534*bbCi^24*dbbi*j2^4*pCi^4*Cos[2*(tt + tti)] + 
         14580*bbCi^27*dpi*j2^4*pCi^5*Cos[2*(tt + tti)] + 
         48600*bbCi^26*dbbi*j2^4*pCi^6*Cos[2*(tt + tti)] - 
         8100*bbCi^29*dpi*j2^4*pCi^7*Cos[2*(tt + tti)] - 
         22275*bbCi^28*dbbi*j2^4*pCi^8*Cos[2*(tt + tti)] + 
         216*bbCi^17*dxi*j2^4*Cos[2*tt + tti] - 2376*bbCi^19*dxi*j2^4*pCi^2*
          Cos[2*tt + tti] + 7560*bbCi^21*dxi*j2^4*pCi^4*Cos[2*tt + tti] - 
         5400*bbCi^23*dxi*j2^4*pCi^6*Cos[2*tt + tti] + 2160*bbCi^16*dbbi*j2^4*
          xCi*Cos[2*tt + tti] - 4752*bbCi^19*dpi*j2^4*pCi*xCi*
          Cos[2*tt + tti] - 28512*bbCi^18*dbbi*j2^4*pCi^2*xCi*
          Cos[2*tt + tti] + 30240*bbCi^21*dpi*j2^4*pCi^3*xCi*
          Cos[2*tt + tti] + 105840*bbCi^20*dbbi*j2^4*pCi^4*xCi*
          Cos[2*tt + tti] - 32400*bbCi^23*dpi*j2^4*pCi^5*xCi*
          Cos[2*tt + tti] - 86400*bbCi^22*dbbi*j2^4*pCi^6*xCi*
          Cos[2*tt + tti] - 432*bbCi^11*j2^3*kC*R*rhoC*yCi*Cos[2*tt + tti] + 
         2592*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*Cos[2*tt + tti] - 
         2160*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*Cos[2*tt + tti] + 
         432*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[2*tt + tti] - 
         2592*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[2*tt + tti] + 
         2160*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[2*tt + tti] + 
         378*bbCi^20*dbbi*j2^4*Cos[3*tt + tti] - 648*bbCi^23*dpi*j2^4*pCi*
          Cos[3*tt + tti] - 5184*bbCi^22*dbbi*j2^4*pCi^2*Cos[3*tt + tti] + 
         4968*bbCi^25*dpi*j2^4*pCi^3*Cos[3*tt + tti] + 22356*bbCi^24*dbbi*
          j2^4*pCi^4*Cos[3*tt + tti] - 9720*bbCi^27*dpi*j2^4*pCi^5*
          Cos[3*tt + tti] - 32400*bbCi^26*dbbi*j2^4*pCi^6*Cos[3*tt + tti] + 
         5400*bbCi^29*dpi*j2^4*pCi^7*Cos[3*tt + tti] + 14850*bbCi^28*dbbi*
          j2^4*pCi^8*Cos[3*tt + tti] - 108*bbCi^17*dxi*j2^4*Cos[tt + 2*tti] + 
         1188*bbCi^19*dxi*j2^4*pCi^2*Cos[tt + 2*tti] - 3780*bbCi^21*dxi*j2^4*
          pCi^4*Cos[tt + 2*tti] + 2700*bbCi^23*dxi*j2^4*pCi^6*
          Cos[tt + 2*tti] - 1080*bbCi^16*dbbi*j2^4*xCi*Cos[tt + 2*tti] + 
         2376*bbCi^19*dpi*j2^4*pCi*xCi*Cos[tt + 2*tti] + 
         14256*bbCi^18*dbbi*j2^4*pCi^2*xCi*Cos[tt + 2*tti] - 
         15120*bbCi^21*dpi*j2^4*pCi^3*xCi*Cos[tt + 2*tti] - 
         52920*bbCi^20*dbbi*j2^4*pCi^4*xCi*Cos[tt + 2*tti] + 
         16200*bbCi^23*dpi*j2^4*pCi^5*xCi*Cos[tt + 2*tti] + 
         43200*bbCi^22*dbbi*j2^4*pCi^6*xCi*Cos[tt + 2*tti] + 
         144*bbCi^11*j2^3*kC*R*rhoC*yCi*Cos[tt + 2*tti] + 
         720*bbCi^13*j2^3*kC*pCi^2*R*rhoC*yCi*Cos[tt + 2*tti] - 
         1440*bbCi^15*j2^3*kC*pCi^4*R*rhoC*yCi*Cos[tt + 2*tti] - 
         144*bbCi^11*j2^3*kD*R*rhoD*yCi*Cos[tt + 2*tti] - 
         720*bbCi^13*j2^3*kD*pCi^2*R*rhoD*yCi*Cos[tt + 2*tti] + 
         1440*bbCi^15*j2^3*kD*pCi^4*R*rhoD*yCi*Cos[tt + 2*tti] + 
         1449*bbCi^20*dbbi*j2^4*Cos[tt + 3*tti] - 2484*bbCi^23*dpi*j2^4*pCi*
          Cos[tt + 3*tti] - 19872*bbCi^22*dbbi*j2^4*pCi^2*Cos[tt + 3*tti] + 
         19044*bbCi^25*dpi*j2^4*pCi^3*Cos[tt + 3*tti] + 85698*bbCi^24*dbbi*
          j2^4*pCi^4*Cos[tt + 3*tti] - 37260*bbCi^27*dpi*j2^4*pCi^5*
          Cos[tt + 3*tti] - 124200*bbCi^26*dbbi*j2^4*pCi^6*Cos[tt + 3*tti] + 
         20700*bbCi^29*dpi*j2^4*pCi^7*Cos[tt + 3*tti] + 56925*bbCi^28*dbbi*
          j2^4*pCi^8*Cos[tt + 3*tti] + 144*bbCi^13*dyi*j2^3*Sin[tt] + 
         864*bbCi^12*dbbi*dyi*j2^3*Sin[tt] + 288*bbCi^17*dyi*j2^4*Sin[tt] - 
         2880*bbCi^15*dpi*dyi*j2^3*pCi*Sin[tt] - 1440*bbCi^15*dyi*j2^3*pCi^2*
          Sin[tt] - 11520*bbCi^14*dbbi*dyi*j2^3*pCi^2*Sin[tt] - 
         3600*bbCi^19*dyi*j2^4*pCi^2*Sin[tt] + 14400*bbCi^17*dpi*dyi*j2^3*
          pCi^3*Sin[tt] + 3600*bbCi^17*dyi*j2^3*pCi^4*Sin[tt] + 
         36000*bbCi^16*dbbi*dyi*j2^3*pCi^4*Sin[tt] + 14400*bbCi^21*dyi*j2^4*
          pCi^4*Sin[tt] - 18000*bbCi^23*dyi*j2^4*pCi^6*Sin[tt] - 
         384*bbCi^7*dxi*j2^2*kD*R*rhoD*Sin[tt] + 1920*bbCi^9*dxi*j2^2*kD*
          pCi^2*R*rhoD*Sin[tt] + 384*bbCi^7*j2^2*kC*R*rhoC*xCi*Sin[tt] - 
         408*bbCi^11*j2^3*kC*R*rhoC*xCi*Sin[tt] - 1920*bbCi^9*j2^2*kC*pCi^2*R*
          rhoC*xCi*Sin[tt] + 3360*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*Sin[tt] - 
         4680*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*Sin[tt] - 
         384*bbCi^7*j2^2*kD*R*rhoD*xCi*Sin[tt] + 408*bbCi^11*j2^3*kD*R*rhoD*
          xCi*Sin[tt] + 3840*bbCi^9*dpi*j2^2*kD*pCi*R*rhoD*xCi*Sin[tt] + 
         1920*bbCi^9*j2^2*kD*pCi^2*R*rhoD*xCi*Sin[tt] + 3840*bbCi^8*dbbi*j2^2*
          kD*pCi^2*R*rhoD*xCi*Sin[tt] - 3360*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*
          Sin[tt] + 4680*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*Sin[tt] + 
         864*bbCi^12*dbbi*j2^3*yCi*Sin[tt] + 2160*bbCi^11*dbbi^2*j2^3*yCi*
          Sin[tt] - 1440*bbCi^15*dpi^2*j2^3*yCi*Sin[tt] + 
         2880*bbCi^16*dbbi*j2^4*yCi*Sin[tt] - 2880*bbCi^15*dpi*j2^3*pCi*yCi*
          Sin[tt] - 23040*bbCi^14*dbbi*dpi*j2^3*pCi*yCi*Sin[tt] - 
         7200*bbCi^19*dpi*j2^4*pCi*yCi*Sin[tt] - 11520*bbCi^14*dbbi*j2^3*
          pCi^2*yCi*Sin[tt] - 40320*bbCi^13*dbbi^2*j2^3*pCi^2*yCi*Sin[tt] + 
         21600*bbCi^17*dpi^2*j2^3*pCi^2*yCi*Sin[tt] - 43200*bbCi^18*dbbi*j2^4*
          pCi^2*yCi*Sin[tt] + 14400*bbCi^17*dpi*j2^3*pCi^3*yCi*Sin[tt] + 
         144000*bbCi^16*dbbi*dpi*j2^3*pCi^3*yCi*Sin[tt] + 
         57600*bbCi^21*dpi*j2^4*pCi^3*yCi*Sin[tt] + 36000*bbCi^16*dbbi*j2^3*
          pCi^4*yCi*Sin[tt] + 162000*bbCi^15*dbbi^2*j2^3*pCi^4*yCi*Sin[tt] + 
         201600*bbCi^20*dbbi*j2^4*pCi^4*yCi*Sin[tt] - 108000*bbCi^23*dpi*j2^4*
          pCi^5*yCi*Sin[tt] - 288000*bbCi^22*dbbi*j2^4*pCi^6*yCi*Sin[tt] + 
         640*bbCi*j2*kC^2*R^2*rhoC^2*yCi*Sin[tt] - 640*bbCi*j2*kD^2*R^2*
          rhoD^2*yCi*Sin[tt] + 486*bbCi^15*j2^3*kC*R*rhoC*Sin[2*tt] - 
         4698*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[2*tt] + 13122*bbCi^19*j2^3*kC*
          pCi^4*R*rhoC*Sin[2*tt] - 8910*bbCi^21*j2^3*kC*pCi^6*R*rhoC*
          Sin[2*tt] - 486*bbCi^15*j2^3*kD*R*rhoD*Sin[2*tt] + 
         4698*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[2*tt] - 13122*bbCi^19*j2^3*kD*
          pCi^4*R*rhoD*Sin[2*tt] + 8910*bbCi^21*j2^3*kD*pCi^6*R*rhoD*
          Sin[2*tt] - 576*bbCi^13*dyi*j2^4*xCi*Sin[2*tt] + 
         5760*bbCi^15*dyi*j2^4*pCi^2*xCi*Sin[2*tt] - 14400*bbCi^17*dyi*j2^4*
          pCi^4*xCi*Sin[2*tt] - 576*bbCi^7*j2^3*kC*R*rhoC*xCi^2*Sin[2*tt] + 
         2880*bbCi^9*j2^3*kC*pCi^2*R*rhoC*xCi^2*Sin[2*tt] + 
         576*bbCi^7*j2^3*kD*R*rhoD*xCi^2*Sin[2*tt] - 2880*bbCi^9*j2^3*kD*
          pCi^2*R*rhoD*xCi^2*Sin[2*tt] - 576*bbCi^13*dxi*j2^4*yCi*Sin[2*tt] + 
         5760*bbCi^15*dxi*j2^4*pCi^2*yCi*Sin[2*tt] - 14400*bbCi^17*dxi*j2^4*
          pCi^4*yCi*Sin[2*tt] - 3456*bbCi^12*dbbi*j2^4*xCi*yCi*Sin[2*tt] + 
         11520*bbCi^15*dpi*j2^4*pCi*xCi*yCi*Sin[2*tt] + 46080*bbCi^14*dbbi*
          j2^4*pCi^2*xCi*yCi*Sin[2*tt] - 57600*bbCi^17*dpi*j2^4*pCi^3*xCi*yCi*
          Sin[2*tt] - 144000*bbCi^16*dbbi*j2^4*pCi^4*xCi*yCi*Sin[2*tt] + 
         576*bbCi^7*j2^3*kC*R*rhoC*yCi^2*Sin[2*tt] - 2880*bbCi^9*j2^3*kC*
          pCi^2*R*rhoC*yCi^2*Sin[2*tt] - 576*bbCi^7*j2^3*kD*R*rhoD*yCi^2*
          Sin[2*tt] + 2880*bbCi^9*j2^3*kD*pCi^2*R*rhoD*yCi^2*Sin[2*tt] - 
         72*bbCi^17*dyi*j2^4*Sin[3*tt] + 792*bbCi^19*dyi*j2^4*pCi^2*
          Sin[3*tt] - 2520*bbCi^21*dyi*j2^4*pCi^4*Sin[3*tt] + 
         1800*bbCi^23*dyi*j2^4*pCi^6*Sin[3*tt] + 264*bbCi^11*j2^3*kC*R*rhoC*
          xCi*Sin[3*tt] - 1440*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*Sin[3*tt] + 
         600*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*Sin[3*tt] - 
         264*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[3*tt] + 1440*bbCi^13*j2^3*kD*
          pCi^2*R*rhoD*xCi*Sin[3*tt] - 600*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*
          Sin[3*tt] - 720*bbCi^16*dbbi*j2^4*yCi*Sin[3*tt] + 
         1584*bbCi^19*dpi*j2^4*pCi*yCi*Sin[3*tt] + 9504*bbCi^18*dbbi*j2^4*
          pCi^2*yCi*Sin[3*tt] - 10080*bbCi^21*dpi*j2^4*pCi^3*yCi*Sin[3*tt] - 
         35280*bbCi^20*dbbi*j2^4*pCi^4*yCi*Sin[3*tt] + 10800*bbCi^23*dpi*j2^4*
          pCi^5*yCi*Sin[3*tt] + 28800*bbCi^22*dbbi*j2^4*pCi^6*yCi*Sin[3*tt] + 
         117*bbCi^15*j2^3*kC*R*rhoC*Sin[tt - 5*tti] + 315*bbCi^17*j2^3*kC*
          pCi^2*R*rhoC*Sin[tt - 5*tti] - 477*bbCi^19*j2^3*kC*pCi^4*R*rhoC*
          Sin[tt - 5*tti] + 45*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[tt - 5*tti] - 
         117*bbCi^15*j2^3*kD*R*rhoD*Sin[tt - 5*tti] - 315*bbCi^17*j2^3*kD*
          pCi^2*R*rhoD*Sin[tt - 5*tti] + 477*bbCi^19*j2^3*kD*pCi^4*R*rhoD*
          Sin[tt - 5*tti] - 45*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[tt - 5*tti] - 
         216*bbCi^11*j2^3*kC*R*rhoC*xCi*Sin[tt - 4*tti] + 
         1296*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*Sin[tt - 4*tti] - 
         1080*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*Sin[tt - 4*tti] + 
         216*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[tt - 4*tti] - 
         1296*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*Sin[tt - 4*tti] + 
         1080*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*Sin[tt - 4*tti] - 
         336*bbCi^11*j2^2*kC*R*rhoC*Sin[tt - 3*tti] + 285*bbCi^15*j2^3*kC*R*
          rhoC*Sin[tt - 3*tti] + 2016*bbCi^13*j2^2*kC*pCi^2*R*rhoC*
          Sin[tt - 3*tti] - 3843*bbCi^17*j2^3*kC*pCi^2*R*rhoC*
          Sin[tt - 3*tti] - 1680*bbCi^15*j2^2*kC*pCi^4*R*rhoC*
          Sin[tt - 3*tti] + 11847*bbCi^19*j2^3*kC*pCi^4*R*rhoC*
          Sin[tt - 3*tti] - 7425*bbCi^21*j2^3*kC*pCi^6*R*rhoC*
          Sin[tt - 3*tti] + 336*bbCi^11*j2^2*kD*R*rhoD*Sin[tt - 3*tti] + 
         1344*bbCi^10*dbbi*j2^2*kD*R*rhoD*Sin[tt - 3*tti] - 
         285*bbCi^15*j2^3*kD*R*rhoD*Sin[tt - 3*tti] - 4032*bbCi^13*dpi*j2^2*
          kD*pCi*R*rhoD*Sin[tt - 3*tti] - 2016*bbCi^13*j2^2*kD*pCi^2*R*rhoD*
          Sin[tt - 3*tti] - 12096*bbCi^12*dbbi*j2^2*kD*pCi^2*R*rhoD*
          Sin[tt - 3*tti] + 3843*bbCi^17*j2^3*kD*pCi^2*R*rhoD*
          Sin[tt - 3*tti] + 6720*bbCi^15*dpi*j2^2*kD*pCi^3*R*rhoD*
          Sin[tt - 3*tti] + 1680*bbCi^15*j2^2*kD*pCi^4*R*rhoD*
          Sin[tt - 3*tti] + 13440*bbCi^14*dbbi*j2^2*kD*pCi^4*R*rhoD*
          Sin[tt - 3*tti] - 11847*bbCi^19*j2^3*kD*pCi^4*R*rhoD*
          Sin[tt - 3*tti] + 7425*bbCi^21*j2^3*kD*pCi^6*R*rhoD*
          Sin[tt - 3*tti] - 441*bbCi^15*j2^3*kC*R*rhoC*Sin[2*(tt - 3*tti)] + 
         3087*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[2*(tt - 3*tti)] - 
         4851*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[2*(tt - 3*tti)] + 
         2205*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[2*(tt - 3*tti)] + 
         441*bbCi^15*j2^3*kD*R*rhoD*Sin[2*(tt - 3*tti)] - 
         3087*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[2*(tt - 3*tti)] + 
         4851*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[2*(tt - 3*tti)] - 
         2205*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[2*(tt - 3*tti)] + 
         504*bbCi^17*dyi*j2^4*Sin[2*tt - 3*tti] - 5544*bbCi^19*dyi*j2^4*pCi^2*
          Sin[2*tt - 3*tti] + 17640*bbCi^21*dyi*j2^4*pCi^4*
          Sin[2*tt - 3*tti] - 12600*bbCi^23*dyi*j2^4*pCi^6*
          Sin[2*tt - 3*tti] + 1008*bbCi^11*j2^3*kC*R*rhoC*xCi*
          Sin[2*tt - 3*tti] - 6048*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*
          Sin[2*tt - 3*tti] + 5040*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*
          Sin[2*tt - 3*tti] - 1008*bbCi^11*j2^3*kD*R*rhoD*xCi*
          Sin[2*tt - 3*tti] + 6048*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*
          Sin[2*tt - 3*tti] - 5040*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*
          Sin[2*tt - 3*tti] + 5040*bbCi^16*dbbi*j2^4*yCi*Sin[2*tt - 3*tti] - 
         11088*bbCi^19*dpi*j2^4*pCi*yCi*Sin[2*tt - 3*tti] - 
         66528*bbCi^18*dbbi*j2^4*pCi^2*yCi*Sin[2*tt - 3*tti] + 
         70560*bbCi^21*dpi*j2^4*pCi^3*yCi*Sin[2*tt - 3*tti] + 
         246960*bbCi^20*dbbi*j2^4*pCi^4*yCi*Sin[2*tt - 3*tti] - 
         75600*bbCi^23*dpi*j2^4*pCi^5*yCi*Sin[2*tt - 3*tti] - 
         201600*bbCi^22*dbbi*j2^4*pCi^6*yCi*Sin[2*tt - 3*tti] - 
         108*bbCi^17*dyi*j2^4*Sin[tt - 2*tti] + 1188*bbCi^19*dyi*j2^4*pCi^2*
          Sin[tt - 2*tti] - 3780*bbCi^21*dyi*j2^4*pCi^4*Sin[tt - 2*tti] + 
         2700*bbCi^23*dyi*j2^4*pCi^6*Sin[tt - 2*tti] - 432*bbCi^11*j2^3*kC*R*
          rhoC*xCi*Sin[tt - 2*tti] + 1008*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*
          Sin[tt - 2*tti] + 432*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[tt - 2*tti] - 
         1008*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*Sin[tt - 2*tti] - 
         1080*bbCi^16*dbbi*j2^4*yCi*Sin[tt - 2*tti] + 2376*bbCi^19*dpi*j2^4*
          pCi*yCi*Sin[tt - 2*tti] + 14256*bbCi^18*dbbi*j2^4*pCi^2*yCi*
          Sin[tt - 2*tti] - 15120*bbCi^21*dpi*j2^4*pCi^3*yCi*
          Sin[tt - 2*tti] - 52920*bbCi^20*dbbi*j2^4*pCi^4*yCi*
          Sin[tt - 2*tti] + 16200*bbCi^23*dpi*j2^4*pCi^5*yCi*
          Sin[tt - 2*tti] + 43200*bbCi^22*dbbi*j2^4*pCi^6*yCi*
          Sin[tt - 2*tti] + 1134*bbCi^15*j2^3*kC*R*rhoC*Sin[2*(tt - 2*tti)] - 
         10962*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[2*(tt - 2*tti)] + 
         30618*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[2*(tt - 2*tti)] - 
         20790*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[2*(tt - 2*tti)] - 
         1134*bbCi^15*j2^3*kD*R*rhoD*Sin[2*(tt - 2*tti)] + 
         10962*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[2*(tt - 2*tti)] - 
         30618*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[2*(tt - 2*tti)] + 
         20790*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[2*(tt - 2*tti)] + 
         288*bbCi^11*j2^2*kC*R*rhoC*Sin[tt - tti] - 27*bbCi^15*j2^3*kC*R*rhoC*
          Sin[tt - tti] - 2304*bbCi^13*j2^2*kC*pCi^2*R*rhoC*Sin[tt - tti] + 
         2799*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[tt - tti] + 
         4320*bbCi^15*j2^2*kC*pCi^4*R*rhoC*Sin[tt - tti] - 
         11205*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[tt - tti] + 
         15345*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[tt - tti] - 
         288*bbCi^11*j2^2*kD*R*rhoD*Sin[tt - tti] - 1152*bbCi^10*dbbi*j2^2*kD*
          R*rhoD*Sin[tt - tti] + 27*bbCi^15*j2^3*kD*R*rhoD*Sin[tt - tti] + 
         4608*bbCi^13*dpi*j2^2*kD*pCi*R*rhoD*Sin[tt - tti] + 
         2304*bbCi^13*j2^2*kD*pCi^2*R*rhoD*Sin[tt - tti] + 
         13824*bbCi^12*dbbi*j2^2*kD*pCi^2*R*rhoD*Sin[tt - tti] - 
         2799*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[tt - tti] - 
         17280*bbCi^15*dpi*j2^2*kD*pCi^3*R*rhoD*Sin[tt - tti] - 
         4320*bbCi^15*j2^2*kD*pCi^4*R*rhoD*Sin[tt - tti] - 
         34560*bbCi^14*dbbi*j2^2*kD*pCi^4*R*rhoD*Sin[tt - tti] + 
         11205*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[tt - tti] - 
         15345*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[tt - tti] - 
         1026*bbCi^15*j2^3*kC*R*rhoC*Sin[2*(tt - tti)] + 
         10422*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[2*(tt - tti)] - 
         34614*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[2*(tt - tti)] + 
         40770*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[2*(tt - tti)] + 
         1026*bbCi^15*j2^3*kD*R*rhoD*Sin[2*(tt - tti)] - 
         10422*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[2*(tt - tti)] + 
         34614*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[2*(tt - tti)] - 
         40770*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[2*(tt - tti)] - 
         231*bbCi^15*j2^3*kC*R*rhoC*Sin[3*(tt - tti)] + 1491*bbCi^17*j2^3*kC*
          pCi^2*R*rhoC*Sin[3*(tt - tti)] - 1785*bbCi^19*j2^3*kC*pCi^4*R*rhoC*
          Sin[3*(tt - tti)] + 525*bbCi^21*j2^3*kC*pCi^6*R*rhoC*
          Sin[3*(tt - tti)] + 231*bbCi^15*j2^3*kD*R*rhoD*Sin[3*(tt - tti)] - 
         1491*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[3*(tt - tti)] + 
         1785*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[3*(tt - tti)] - 
         525*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[3*(tt - tti)] - 
         432*bbCi^17*dyi*j2^4*Sin[2*tt - tti] + 5616*bbCi^19*dyi*j2^4*pCi^2*
          Sin[2*tt - tti] - 23760*bbCi^21*dyi*j2^4*pCi^4*Sin[2*tt - tti] + 
         32400*bbCi^23*dyi*j2^4*pCi^6*Sin[2*tt - tti] - 1296*bbCi^11*j2^3*kC*
          R*rhoC*xCi*Sin[2*tt - tti] + 11232*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*
          Sin[2*tt - tti] - 23760*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*
          Sin[2*tt - tti] + 1296*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[2*tt - tti] - 
         11232*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*Sin[2*tt - tti] + 
         23760*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*Sin[2*tt - tti] - 
         4320*bbCi^16*dbbi*j2^4*yCi*Sin[2*tt - tti] + 11232*bbCi^19*dpi*j2^4*
          pCi*yCi*Sin[2*tt - tti] + 67392*bbCi^18*dbbi*j2^4*pCi^2*yCi*
          Sin[2*tt - tti] - 95040*bbCi^21*dpi*j2^4*pCi^3*yCi*
          Sin[2*tt - tti] - 332640*bbCi^20*dbbi*j2^4*pCi^4*yCi*
          Sin[2*tt - tti] + 194400*bbCi^23*dpi*j2^4*pCi^5*yCi*
          Sin[2*tt - tti] + 518400*bbCi^22*dbbi*j2^4*pCi^6*yCi*
          Sin[2*tt - tti] + 198*bbCi^15*j2^3*kC*R*rhoC*Sin[3*tt - tti] - 
         1674*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[3*tt - tti] + 
         3690*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[3*tt - tti] - 
         1350*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[3*tt - tti] - 
         198*bbCi^15*j2^3*kD*R*rhoD*Sin[3*tt - tti] + 1674*bbCi^17*j2^3*kD*
          pCi^2*R*rhoD*Sin[3*tt - tti] - 3690*bbCi^19*j2^3*kD*pCi^4*R*rhoD*
          Sin[3*tt - tti] + 1350*bbCi^21*j2^3*kD*pCi^6*R*rhoD*
          Sin[3*tt - tti] - 144*bbCi^11*j2^3*kC*R*rhoC*xCi*Sin[tti] + 
         1440*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*Sin[tti] - 
         3600*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*Sin[tti] + 
         144*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[tti] - 1440*bbCi^13*j2^3*kD*pCi^2*
          R*rhoD*xCi*Sin[tti] + 3600*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*
          Sin[tti] - 72*bbCi^15*j2^3*kC*R*rhoC*Sin[2*tti] + 
         792*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[2*tti] - 2520*bbCi^19*j2^3*kC*
          pCi^4*R*rhoC*Sin[2*tti] + 1800*bbCi^21*j2^3*kC*pCi^6*R*rhoC*
          Sin[2*tti] + 72*bbCi^15*j2^3*kD*R*rhoD*Sin[2*tti] - 
         792*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[2*tti] + 2520*bbCi^19*j2^3*kD*
          pCi^4*R*rhoD*Sin[2*tti] - 1800*bbCi^21*j2^3*kD*pCi^6*R*rhoD*
          Sin[2*tti] - 144*bbCi^11*j2^2*kC*R*rhoC*Sin[tt + tti] + 
         207*bbCi^15*j2^3*kC*R*rhoC*Sin[tt + tti] + 864*bbCi^13*j2^2*kC*pCi^2*
          R*rhoC*Sin[tt + tti] - 3537*bbCi^17*j2^3*kC*pCi^2*R*rhoC*
          Sin[tt + tti] - 720*bbCi^15*j2^2*kC*pCi^4*R*rhoC*Sin[tt + tti] + 
         13365*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[tt + tti] - 
         10035*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[tt + tti] + 
         144*bbCi^11*j2^2*kD*R*rhoD*Sin[tt + tti] + 576*bbCi^10*dbbi*j2^2*kD*
          R*rhoD*Sin[tt + tti] - 207*bbCi^15*j2^3*kD*R*rhoD*Sin[tt + tti] - 
         1728*bbCi^13*dpi*j2^2*kD*pCi*R*rhoD*Sin[tt + tti] - 
         864*bbCi^13*j2^2*kD*pCi^2*R*rhoD*Sin[tt + tti] - 
         5184*bbCi^12*dbbi*j2^2*kD*pCi^2*R*rhoD*Sin[tt + tti] + 
         3537*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[tt + tti] + 
         2880*bbCi^15*dpi*j2^2*kD*pCi^3*R*rhoD*Sin[tt + tti] + 
         720*bbCi^15*j2^2*kD*pCi^4*R*rhoD*Sin[tt + tti] + 
         5760*bbCi^14*dbbi*j2^2*kD*pCi^4*R*rhoD*Sin[tt + tti] - 
         13365*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[tt + tti] + 
         10035*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[tt + tti] - 
         81*bbCi^15*j2^3*kC*R*rhoC*Sin[2*(tt + tti)] + 567*bbCi^17*j2^3*kC*
          pCi^2*R*rhoC*Sin[2*(tt + tti)] - 891*bbCi^19*j2^3*kC*pCi^4*R*rhoC*
          Sin[2*(tt + tti)] + 405*bbCi^21*j2^3*kC*pCi^6*R*rhoC*
          Sin[2*(tt + tti)] + 81*bbCi^15*j2^3*kD*R*rhoD*Sin[2*(tt + tti)] - 
         567*bbCi^17*j2^3*kD*pCi^2*R*rhoD*Sin[2*(tt + tti)] + 
         891*bbCi^19*j2^3*kD*pCi^4*R*rhoD*Sin[2*(tt + tti)] - 
         405*bbCi^21*j2^3*kD*pCi^6*R*rhoD*Sin[2*(tt + tti)] + 
         216*bbCi^17*dyi*j2^4*Sin[2*tt + tti] - 2376*bbCi^19*dyi*j2^4*pCi^2*
          Sin[2*tt + tti] + 7560*bbCi^21*dyi*j2^4*pCi^4*Sin[2*tt + tti] - 
         5400*bbCi^23*dyi*j2^4*pCi^6*Sin[2*tt + tti] + 432*bbCi^11*j2^3*kC*R*
          rhoC*xCi*Sin[2*tt + tti] - 2592*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*
          Sin[2*tt + tti] + 2160*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*
          Sin[2*tt + tti] - 432*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[2*tt + tti] + 
         2592*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*Sin[2*tt + tti] - 
         2160*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*Sin[2*tt + tti] + 
         2160*bbCi^16*dbbi*j2^4*yCi*Sin[2*tt + tti] - 4752*bbCi^19*dpi*j2^4*
          pCi*yCi*Sin[2*tt + tti] - 28512*bbCi^18*dbbi*j2^4*pCi^2*yCi*
          Sin[2*tt + tti] + 30240*bbCi^21*dpi*j2^4*pCi^3*yCi*
          Sin[2*tt + tti] + 105840*bbCi^20*dbbi*j2^4*pCi^4*yCi*
          Sin[2*tt + tti] - 32400*bbCi^23*dpi*j2^4*pCi^5*yCi*
          Sin[2*tt + tti] - 86400*bbCi^22*dbbi*j2^4*pCi^6*yCi*
          Sin[2*tt + tti] - 99*bbCi^15*j2^3*kC*R*rhoC*Sin[3*tt + tti] + 
         639*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[3*tt + tti] - 
         765*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[3*tt + tti] + 
         225*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[3*tt + tti] + 
         99*bbCi^15*j2^3*kD*R*rhoD*Sin[3*tt + tti] - 639*bbCi^17*j2^3*kD*
          pCi^2*R*rhoD*Sin[3*tt + tti] + 765*bbCi^19*j2^3*kD*pCi^4*R*rhoD*
          Sin[3*tt + tti] - 225*bbCi^21*j2^3*kD*pCi^6*R*rhoD*
          Sin[3*tt + tti] - 108*bbCi^17*dyi*j2^4*Sin[tt + 2*tti] + 
         1188*bbCi^19*dyi*j2^4*pCi^2*Sin[tt + 2*tti] - 3780*bbCi^21*dyi*j2^4*
          pCi^4*Sin[tt + 2*tti] + 2700*bbCi^23*dyi*j2^4*pCi^6*
          Sin[tt + 2*tti] - 144*bbCi^11*j2^3*kC*R*rhoC*xCi*Sin[tt + 2*tti] - 
         720*bbCi^13*j2^3*kC*pCi^2*R*rhoC*xCi*Sin[tt + 2*tti] + 
         1440*bbCi^15*j2^3*kC*pCi^4*R*rhoC*xCi*Sin[tt + 2*tti] + 
         144*bbCi^11*j2^3*kD*R*rhoD*xCi*Sin[tt + 2*tti] + 
         720*bbCi^13*j2^3*kD*pCi^2*R*rhoD*xCi*Sin[tt + 2*tti] - 
         1440*bbCi^15*j2^3*kD*pCi^4*R*rhoD*xCi*Sin[tt + 2*tti] - 
         1080*bbCi^16*dbbi*j2^4*yCi*Sin[tt + 2*tti] + 2376*bbCi^19*dpi*j2^4*
          pCi*yCi*Sin[tt + 2*tti] + 14256*bbCi^18*dbbi*j2^4*pCi^2*yCi*
          Sin[tt + 2*tti] - 15120*bbCi^21*dpi*j2^4*pCi^3*yCi*
          Sin[tt + 2*tti] - 52920*bbCi^20*dbbi*j2^4*pCi^4*yCi*
          Sin[tt + 2*tti] + 16200*bbCi^23*dpi*j2^4*pCi^5*yCi*
          Sin[tt + 2*tti] + 43200*bbCi^22*dbbi*j2^4*pCi^6*yCi*
          Sin[tt + 2*tti] + 210*bbCi^15*j2^3*kC*R*rhoC*Sin[tt + 3*tti] - 
         1422*bbCi^17*j2^3*kC*pCi^2*R*rhoC*Sin[tt + 3*tti] + 
         2502*bbCi^19*j2^3*kC*pCi^4*R*rhoC*Sin[tt + 3*tti] - 
         1290*bbCi^21*j2^3*kC*pCi^6*R*rhoC*Sin[tt + 3*tti] - 
         210*bbCi^15*j2^3*kD*R*rhoD*Sin[tt + 3*tti] + 1422*bbCi^17*j2^3*kD*
          pCi^2*R*rhoD*Sin[tt + 3*tti] - 2502*bbCi^19*j2^3*kD*pCi^4*R*rhoD*
          Sin[tt + 3*tti] + 1290*bbCi^21*j2^3*kD*pCi^6*R*rhoD*
          Sin[tt + 3*tti])))/(196608*bbCi^7)
