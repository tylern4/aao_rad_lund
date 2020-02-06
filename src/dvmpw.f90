!*******************Valery Model****************************************
subroutine DVMPW(COST, W, Q2, Phi_g, E, heli, MESONMASS, &
        sigma0, S_T, S_L, S_TT, S_LT, S_LTP)
    !
    !  dsigma/dQ2 dphi dCosTheta* dW for ep-->ep pi0
    !
    !  exc pi0 x-section
    !
    !input:
    ! COST = Cos(Theta*), Theta* is the angle (pi0-gamma* or  p-p')in the CM system
    ! W is W
    ! xb,Q2 x and Q^2 (GeV^2)
    ! Phi_g angle in the photon frame (radians)
    ! E energy of the electron in GeV
    ! heli electron helicity -1.0 or +1.0
    !
    !del2=t (negative GeV^2)           NEGATIVE !!!!
    !
    ! MESONMASS is the mass of the pi0 or eta.
    ! The actual masses that will be used for the calculations are in pi0eta.par file

    IMPLICIT NONE
    REAL    COST, W, del2, xb, Q2, Phi_g, E, heli, mesonmass
    REAL    Mp, mele, pi
    parameter (Mp = 0.93827)
    parameter (mele = 0.000511)
    parameter (pi = 3.1415926536)

    REAL     sigma0, S_T, S_L, S_LT, S_TT, S_LTP
    REAL     SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT, SIGMA_LTP
    REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
    EXTERNAL EPSILON, SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT
    EXTERNAL SIGMA_LTP
    REAL     JACG, JACR
    LOGICAL  CHECK_KINE
    !      REAL DVMPX
    INCLUDE 'pi0eta_new.par'
    !

    CALL XSINIT(MESONMASS)

    !
    S_T = 0.0
    S_L = 0.0
    S_LT = 0.0
    S_TT = 0.0
    S_LTP = 0.0
    IF(.NOT.CHECK_KINE(COST, W, Q2, E, del2, xb, JACG, JACR)) RETURN
    EPS = EPSILON(XB, q2, e)
    S_T = 0.001 / (2. * PI) * SIGMA_T   (COST, W, Q2, E)
    S_L = 0.001 / (2. * PI) * SIGMA_L   (COST, W, Q2, E)
    S_LT = 0.001 / (2. * PI) * 2.0 * SIGMA_LT  (COST, W, Q2, E)
    S_TT = 0.001 / (2. * PI) * SIGMA_TT  (COST, W, Q2, E)
    S_LTP = 0.001 / (2. * PI) * 2.0 * SIGMA_LTP (COST, W, Q2, E)

    !
    !       sigma0 =0.000001*FLUXW(xb,Q2,E)/(2.*PI)*JACG*(
    sigma0 = 1 / (2. * PI) * (&
            S_T + EPS * S_L + &
                    EPS * S_TT * COS(2 * PHI_G) + &
                    SQRT(2. * EPS * (1. + EPS)) * S_LT * COS(PHI_G) + &
                    HELI * SQRT(2. * EPS * (1. - EPS)) * S_LTP * SIN(2 * PHI_G)&
            )

    IF (sigma0.LT.0.0) then
        sigma0 = 0.0
        S_T = 0.0
        S_L = 0.0
        S_LT = 0.0
        S_TT = 0.0
        S_LTP = 0.0
    ENDIF
    !       PRINT *,DVMPW,COST,W,Q2,del2,xb,phi_g,JACG,JACR
    !      PRINT *,DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)*JACG*JACR
    !      PRINT *,DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
    !      PRINT *,del2,xb,Q2, Phi_g,E,heli,MESONMASS

    !       PRINT *,DVMPW,COST,W,Q2,phi_g,eps
    !       print *,'DVMP',sigma0, S_T,S_L,S_TT,S_LT
    !       print *,'S_TT  ',  S_TT, EPS     *S_TT  * COS(2*PHI_G)
    !       print *, 'S_LT ',S_LT,SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)

    RETURN
END

!====================================================================================c

REAL FUNCTION SIGMA_T(COST, W, Q2, E)
    IMPLICIT NONE
    REAL          COST, W, del2, x, Q2, E
    EXTERNAL      tminq
    REAL          tminq, T0, SLOPE, T, JACG, JACR
    LOGICAL       CHECK_KINE
    EXTERNAL      XSIGMA_T
    REAL          XSIGMA_T

    IF(CHECK_KINE(COST, W, Q2, E, del2, x, JACG, JACR)) THEN
        SIGMA_T = JACR * XSIGMA_T(del2, x, Q2, E)
    ELSE
        SIGMA_T = 0.0
    ENDIF
    RETURN
END
!=============================================================      
REAL FUNCTION SIGMA_L(COST, W, Q2, E)
    IMPLICIT NONE
    REAL          COST, W, del2, x, Q2, E
    EXTERNAL      tminq
    REAL          tminq, T0, SLOPE, T, JACG, JACR
    LOGICAL CHECK_KINE
    EXTERNAL      XSIGMA_L
    REAL          XSIGMA_L

    IF(CHECK_KINE(COST, W, Q2, E, del2, x, JACG, JACR)) THEN
        SIGMA_L = JACR * XSIGMA_L(del2, x, Q2, E)
    ELSE
        SIGMA_L = 0.0
    ENDIF
    RETURN
END

!=============================================================      
REAL FUNCTION SIGMA_TT(COST, W, Q2, E)
    IMPLICIT NONE
    REAL          COST, W, del2, x, Q2, E
    EXTERNAL      tminq
    REAL          tminq, T0, SLOPE, T, JACG, JACR
    LOGICAL CHECK_KINE
    EXTERNAL      XSIGMA_TT
    REAL          XSIGMA_TT

    IF(CHECK_KINE(COST, W, Q2, E, del2, x, JACG, JACR)) THEN
        SIGMA_TT = JACR * XSIGMA_TT(del2, x, Q2, E)
    ELSE
        SIGMA_TT = 0.0
    ENDIF
    RETURN
END

!=============================================================      
REAL FUNCTION SIGMA_LT(COST, W, Q2, E)
    IMPLICIT NONE
    REAL          COST, W, del2, x, Q2, E
    EXTERNAL      tminq
    REAL          tminq, T0, SLOPE, T, JACG, JACR
    LOGICAL CHECK_KINE
    EXTERNAL      XSIGMA_LT
    REAL          XSIGMA_LT

    IF(CHECK_KINE(COST, W, Q2, E, del2, x, JACG, JACR)) THEN
        SIGMA_LT = JACR * XSIGMA_LT(del2, x, Q2, E)
    ELSE
        SIGMA_LT = 0.0
    ENDIF
    RETURN
END

!====================================================================================c

REAL FUNCTION SIGMA_LTP(COST, W, Q2, E)
    IMPLICIT NONE
    REAL COST, W, del2, x, Q2, E
    REAL EPS, EPSILON, A2_PHI, FLUXW, SIGMA_TOT
    EXTERNAL EPSILON, FLUXW

    SIGMA_LTP = 0.

    RETURN
END

!====================================================================================c

LOGICAL FUNCTION CHECK_KINE(COST, W, Q2, E, del2, xb, JACG, JACR)
    IMPLICIT NONE
    REAL COST, W, del2, xb, Q2, E, JACG, JACR
    real Mp, mele, pi, mpi0
    real fluxw, rxs
    parameter (Mp = 0.93827)
    parameter (mele = 0.000511)
    parameter (pi = 3.1415926536)

    real nu, W2, qmod, E1cm, P1cm, E2cm, P2cm, del2max, del2min
    real  xmin1, xmax1
    real y, e1, epsilon
    INCLUDE 'pi0eta_new.par'

    MPI0 = AM(K)

    CHECK_KINE = .FALSE.
    W2 = W * W
    XB = Q2 / (W2 + Q2 - Mp * Mp)
    nu = (W2 + Q2 - Mp * Mp) / (2 * Mp)
    qmod = sqrt(nu**2 + Q2)
    IF(W         .LT. Mp + Mpi0)              RETURN
    IF(W         .GT. -Q2 + 2 * Mp * E + Mp * Mp)     RETURN
    IF(Q2 / (4 * E * (E - NU)) .GE. 1.0)            RETURN
    IF(XB        .LT. Q2 / (2 * Mp * E))          RETURN
    IF(ABS(COST) .GT. 1.0)                  RETURN
    xmin1 = Q2 / (2.0 * Mp * E)
    xmax1 = 1.0

    E1cm = Mp * (Mp + nu) / W
    E2cm = (W2 + Mp**2 - Mpi0**2) / (2D0 * W)
    IF(E1cm.LE.Mp .OR. E2cm.LE. Mp) RETURN
    P1cm = Mp * qmod / W
    P2cm = SQRT(E2CM**2 - Mp**2)
    del2max = 2.0 * (Mp**2 - E1cm * E2cm - P1cm * P2cm)
    del2min = 2.0 * (Mp**2 - E1cm * E2cm + P1cm * P2cm)
    del2 = del2min - 2. * P1cm * P2cm * (1 - COST)
    IF(xb.le.xmin1 .or. xb.gt.xmax1)         return   !    x  out of range
    IF(del2.ge.del2min .or. del2.le.del2max) return   ! delta out of range

    y = Q2 / (2 * Mp * xb * E)
    e1 = (y * xb * Mp)**2 / Q2
    EPSILON = (1.0 - y - e1) / (1 - y + y**2 / 2 + e1)

    IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

    !      jacobian
    !      JAC=4*P1cm*P2CM*W*XB/(W2+Q2-Mp*Mp)
    JACG = 2 * W * XB / (W2 + Q2 - Mp * Mp)
    JACR = 2 * P1cm * P2CM

    CHECK_KINE = .TRUE.

    RETURN
END

!=======================================================================================C

SUBROUTINE XSINIT(MESONMASS)
    IMPLICIT NONE
    REAL MESONMASS
    INCLUDE 'pi0eta.par'

    IF(MESONMASS.GT.0.140) THEN
        K=2
    ELSE
        K=1
    ENDIF

    RETURN
END


LOGICAL FUNCTION XCHECK_KINE(del2,xb,Q2,E)
    IMPLICIT NONE
    REAL del2,xb,Q2,E
    real Mp, mele, pi,mpi0
    real fluxw, rxs
    parameter (Mp=0.93827)
    parameter (mele=0.000511)
    parameter (pi=3.1415926536)


    real nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
    real  xmin1,xmax1
    real y, e1, epsilon

    INCLUDE 'pi0eta.par'

    MPI0=AM(K)
    XCHECK_KINE=.FALSE.
    xmin1 = Q2/(2.0*Mp*E)
    xmax1 = 1.0
    nu  = Q2/(2D0*Mp*xb)
    W2  = Mp**2 + 2.0*Mp*nu - Q2
    IF(W2.LT.(Mp+Mpi0)**2)                               RETURN
    W   = sqrt(W2)
    qmod = sqrt(nu**2 + Q2)

    E1cm = Mp*(Mp + nu)/W
    P1cm = Mp*qmod/W
    E2cm = (W2 + Mp**2-Mpi0**2)/(2.*W)
    IF(E2cm.LE.Mp)                             RETURN
    P2cm = SQRT(E2CM**2 - Mp**2)
    del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
    del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)

    IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range
    IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

    y=Q2/(2*Mp*xb*E)
    e1=(y*xb*Mp)**2/Q2
    EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)

    IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

    XCHECK_KINE=.TRUE.

    RETURN
END

REAL FUNCTION XSIGMA_LT(del2,x,Q2,E)
    IMPLICIT NONE
    REAL          del2,x,Q2,E
    EXTERNAL      tminq
    REAL          tminq,T0,SLOPE,T
    LOGICAL XCHECK_KINE
    INCLUDE 'pi0eta.par'

    IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        SLOPE = 2.*1.1*ALOG(X)
        T=-DEL2
        XSIGMA_LT=P(9,k)*SQRT(T-T0)*EXP(SLOPE*T*P(10,k))*X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
    ELSE
        XSIGMA_LT=0.0
    ENDIF
    RETURN
END

REAL FUNCTION XSIGMA_LTP(del2,x,Q2,E)
    IMPLICIT NONE
    REAL del2,x,Q2,E
    REAL EPS, EPSILON,A2_PHI,FLUXW,SIGMA_TOT
    EXTERNAL EPSILON,FLUXW
    INCLUDE 'pi0eta.par'

    XSIGMA_LTP=0.

    RETURN
END

REAL FUNCTION tminq(Q2,X)
    IMPLICIT NONE
    REAL Q2,X
    REAL W2,W,S,E1CM,E3CM,P1CM,P3CM,TMAX
    INTEGER I
    real alpha,Mp,PI,Me
    REAL MPI0

    INCLUDE 'pi0eta.par'

    MPI0=AM(K)
    tminq=0.

    Mp=0.9382723
    Me=0.0051
    PI=3.14151926

    IF(X.LE.0. .OR. X.GE.1.)        RETURN
    W2 = Q2*(1./X-1.)+Mp**2
    W=SQRT(W2)
    IF(W.LT.Mp+Mpi0)                RETURN

    E1CM=(W2+Q2+Mp**2)/(2*W)
    P1CM=SQRT(E1CM**2-MP**2)

    E3CM=(W2-MPI0**2+Mp**2)/(2*W)
    P3CM=SQRT(E3CM**2-MP**2)

    TMINQ=-((Q2+MPI0**2)**2/4./W2-(P1CM-P3CM)**2)

    RETURN
END


REAL FUNCTION XSIGMA_L(del2,x,Q2,E)
    IMPLICIT NONE
    REAL          del2,x,Q2,E
    EXTERNAL      tminq
    REAL          tminq,T0,SLOPE,T
    LOGICAL XCHECK_KINE
    INCLUDE 'pi0eta.par'

    IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        SLOPE = 2.*1.1*ALOG(X)
        T=-DEL2
        XSIGMA_L=Q2*(P(4,k)+P(5,k)*(T-T0))*EXP(SLOPE*T*P(6,k))*X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
    ELSE
        XSIGMA_L=0.0
    ENDIF
    RETURN
END

REAL FUNCTION XSIGMA_T(del2,x,Q2,E)
    IMPLICIT NONE
    REAL          del2,x,Q2,E
    EXTERNAL      tminq
    REAL          tminq,T0,SLOPE,T
    LOGICAL       XCHECK_KINE
    INCLUDE 'pi0eta.par'

    IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        SLOPE = 2.*1.1*ALOG(X)
        T=-DEL2
        XSIGMA_T=(P(1,k)+P(2,k)*SQRT(T-T0))*EXP(SLOPE*T*P(3,k))*X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
    ELSE
        XSIGMA_T=0.0
    ENDIF
    RETURN
END

REAL FUNCTION XSIGMA_TT(del2,x,Q2,E)
    IMPLICIT NONE
    REAL          del2,x,Q2,E
    EXTERNAL      tminq
    REAL          tminq,T0,SLOPE,T
    LOGICAL XCHECK_KINE
    INCLUDE 'pi0eta.par'

    IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        SLOPE = 2.*1.1*ALOG(X)
        T=-DEL2
        XSIGMA_TT=P(7,k)*(T-T0)*EXP(SLOPE*T*P(8,k))*X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
    ELSE
        XSIGMA_TT=0.0
    ENDIF
    RETURN
END

real function EPSILON(x,Q2,E)
    implicit none
    real x,Q2,E,y,eps,e1
    real alpha,Mp,PI
    parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)

    y=Q2/(2*Mp*x*E)
    e1=(y*x*Mp)**2/Q2
    EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
    return
end
