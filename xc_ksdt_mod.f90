      MODULE xc_ksdt_mod
!----------------------------------------------------------------------
!
!   DOUBLE PRECISION FUNCTION exc_ksdt(rs,th,xi,exc0,dexc0)
!
!   INPUT: rs = Wigner-Seitz radius, scaled by the Bohr radius
!          t = Ratio of temperature to Fermi temperature
!          xi = Spin-polarization (1 = polarized, 0 = unpolarized)
!
!   RETURNS: exc_ksdt = Exchange-correlation energy (Rydbergs)
!
!
!   DOUBLE PRECISION FUNCTION dexc_ksdt(rs,th,xi,exc0,dexc0)
!
!   INPUT: rs = Wigner-Seitz radius, scaled by the Bohr radius
!          t = Ratio of temperature to Fermi temperature
!          xi = Spin-polarization (1 = polarized, 0 = unpolarized)
!
!   RETURNS: dexc_ksdt = rs derivative of exc_ksdt (Rydbergs/rs)
!
!   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   EXAMPLE USAGE:
!
!      PROGRAM main
!
!         USE xc_ksdt_mod
!         IMPLICIT NONE
!         DOUBLE PRECISION, PARAMETER :: rs = 4.0, t = 1.0, xi = 1
!         DOUBLE PRECISION :: exc, dexc
!
!         exc = exc_ksdt(rs,t,xi)
!         dexc = dexc_ksdt(rs,t,xi)
!         write(*,*) exc, dexc
!
!      END PROGRAM main
!
!   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   If you use this module in your calculations, please cite both the
!   fit and original Monte Carlo simulation:
!
!      V. V. Karasiev, T. Sjostrom, J. Dufty, S.B. Trickey
!      Local Density Approximation Exchange-correlation Free-energy
!      Functional
!      arXiv:1311.4903
!
!      E. W. Brown, B. K. Clark, J. L. DuBois, and D. M. Ceperley
!      Path Integral Monte Carlo simulation of the warm-dense
!      homogeneous electron gas
!      Phys. Rev. Lett. 110, 146405 (2013)
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      CONTAINS

        SUBROUTINE params(xi,w,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,  &
                          c1,c2,c3,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)
          DOUBLE PRECISION :: xi
          DOUBLE PRECISION :: w,a1,a2,a3,a4,a5,a6,a7
          DOUBLE PRECISION :: b1,b2,b3,b4,b5,c1,c2,c3
          DOUBLE PRECISION :: d1,d2,d3,d4,d5,e1,e2,e3,e4,e5
          a1 = 0.610887
          a2 = 0.75
          a3 = 3.04363
          a4 = -0.09227
          a5 = 1.7035
          a6 = 8.31051
          a7 = 5.1105
          IF (xi == 0) THEN
            w = 1.
            b1 = 0.283991
            b2 = 48.934733
            b3 = 0.371003
            b4 = 61.093930
            b5 = sqrt(3.)*b3
            c1 = 0.870091
            c2 = 0.193157
            c3 = 2.415881
            d1 = 0.579823
            d2 = 94.541433
            d3 = 97.834705
            d4 = 59.950878
            d5 = 24.385571
            e1 = 0.212027
            e2 = 16.720565
            e3 = 28.488645
            e4 = 34.032003
            e5 = 17.238061
          ELSE
            w = 2.**(1./3.)
            b1 = 0.32909
            b2 = 86.234018
            b3 = 0.711881
            b4 = 104.908343
            b5 = sqrt(3.)*b3
            c1 = 0.849303
            c2 = 0.172616
            c3 = 0.093613
            d1 = 0.551387
            d2 = 141.454167
            d3 = 91.832049
            d4 = 106.726616
            d5 = 21.222501
            e1 = 0.152664
            e2 = 15.087590
            e3 = 25.270051
            e4 = 116.817139
            e5 = 17.295352
          END IF
        END SUBROUTINE params

        SUBROUTINE coefficients(xi,t,w,a,b,c,d,e,da,db,dc,dd,de)
          DOUBLE PRECISION :: xi,t
          DOUBLE PRECISION :: a1,a2,a3,a4,a5,a6,a7
          DOUBLE PRECISION :: b1,b2,b3,b4,b5,c1,c2,c3
          DOUBLE PRECISION :: d1,d2,d3,d4,d5,e1,e2,e3,e4,e5
          DOUBLE PRECISION :: w,a,b,c,d,e,da,db,dc,dd,de

          call params(xi,w,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,      &
                      c1,c2,c3,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)

          a = a1*tanh(1./t)*(a2 + a3*(t**2.) + a4*(t**3.)            &
              + a5*(t**4.))/(1. + a6*(t**2.) + a7*(t**4.))
          b = tanh(1./sqrt(t))*(b1 + b2*(t**2.) + b3*(t**4.))/(1.    &
              + b4*(t**2.) + b5*(t**4.))
          e = tanh(1./t)*(e1 + e2*(t**2.) + e3*(t**4.))/(1.          &
              + e4*(t**2.) + e5*(t**4.))
          c = (c1 + c2*exp(-c3/t))*e
          d = tanh(1./sqrt(t))*(d1 + d2*(t**2.) + d3*(t**4.))/(1.    &
              + d4*(t**2.) + d5*(t**4.))

          da = -((1./(cosh(1./t)*t)**2.)*(a/tanh(1./t)))             &
               + (a1*tanh(1./t)*(2.*a3*t + 3*a4*(t**2.)              &
               + 4.*a5*(t**3.))/(1. + a6*(t**2.) + a7*(t**4.)))      &
               - (a*(2.*a6*t + 4.*a7*(t**3.))/(1. + a6*(t**2.)       &
               + a7*(t**4.)))
          db = -(((1./cosh(1./sqrt(t))**2.)/(2.*(t**1.5)))           &
               * b/tanh(1./sqrt(t)))                                 &
               + (tanh(1./sqrt(t))*(2.*b2*t + 4.*b3*(t**3.))/(1.     &
               + b4*(t**2.) + b5*(t**4.)))                           &
               - (b*(2.*b4*t + 4.*b5*(t**3.))/(1. + b4*(t**2.)       &
               + b5*(t**4.)))
          de = -(((1./cosh(1./t)**2.)/(t**2.))*e/tanh(1./t))         &
               + (tanh(1./t)*(2.*e2*t + 4.*e3*(t**3.))/(1.           &
               + e4*(t**2.) + e5*(t**4.)))                           &
               - (e*(2.*e4*t + 4.*e5*(t**3.))/(1. + e4*(t**2.)       &
               + e5*(t**4.)))
          dc = (c1 + c2*exp(-c3/t))*de                               &
               + (c2*c3*exp(-c3/t)/(t**2.))*e
          dd = -(((1./cosh(1./sqrt(t))**2.)/(2.*(t**1.5)))           &
               * d/tanh(1./sqrt(t)))                                 &
               + (tanh(1./sqrt(t))*(2.*d2*t + 4.*d3*(t**3.))/(1.     &
               + d4*(t**2.) + d5*(t**4.)))                           &
               - (d*(2.*d4*t + 4.*d5*(t**3.))/(1. + d4*(t**2.)       &
               + d5*(t**4.)))
        END SUBROUTINE coefficients

        DOUBLE PRECISION FUNCTION exc0_ksdt(rs,xi)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs,xi
          DOUBLE PRECISION :: w,a1,a2,a3,a4,a5,a6,a7
          DOUBLE PRECISION :: b1,b2,b3,b4,b5,c1,c2,c3
          DOUBLE PRECISION :: d1,d2,d3,d4,d5,e1,e2,e3,e4,e5
          DOUBLE PRECISION :: dem

          call params(xi,w,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,      &
                      c1,c2,c3,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)

          dem = 1. + d1*sqrt(rs) + e1*rs
          exc0_ksdt = (-1./rs)*(w*a1*a2 + b1*sqrt(rs) + c1*e1*rs)/dem
          exc0_ksdt = 2.*exc0_ksdt
        END FUNCTION exc0_ksdt

        DOUBLE PRECISION FUNCTION exc_ksdt(rs,t,xi)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs,t,xi
          DOUBLE PRECISION :: w,a,b,c,d,e,da,db,dc,dd,de
          DOUBLE PRECISION :: dem,fxc,dfxcdt

          IF (t == 0) THEN
            exc_ksdt = exc0_ksdt(rs,xi)
          ELSE
            call coefficients(xi,t,w,a,b,c,d,e,da,db,dc,dd,de)

            dem = 1.+d*sqrt(rs)+e*rs
            fxc = (-1./rs)*(w*a + b*sqrt(rs) + c*rs)/dem
            dfxcdt = - ((1./rs)*(w*da + db*sqrt(rs) + dc*rs)/dem)    &
                     - (fxc*(dd*sqrt(rs) + de*rs)/dem)
            exc_ksdt = 2.*(fxc - t*dfxcdt)
          END IF
        END FUNCTION exc_ksdt

        DOUBLE PRECISION FUNCTION dexc0_ksdt(rs,xi)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs,xi
          DOUBLE PRECISION :: w,a1,a2,a3,a4,a5,a6,a7
          DOUBLE PRECISION :: b1,b2,b3,b4,b5,c1,c2,c3
          DOUBLE PRECISION :: d1,d2,d3,d4,d5,e1,e2,e3,e4,e5
          DOUBLE PRECISION :: dem, exc0, dexc0drs

          call params(xi,w,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,      &
                      c1,c2,c3,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)

          dem = 1. + d1*sqrt(rs) + e1*rs
          exc0 = (-1./rs)*(w*a1*a2 + b1*sqrt(rs) + c1*e1*rs)/dem
          dexc0drs = - (exc0*((d1/(2.*sqrt(rs))) + e1)/dem)          &
                    - ((1./rs)*exc0)                                 &
                    - ((1./rs)*((b1/(2.*sqrt(rs))) + c1*e1)/dem)
          dexc0_ksdt = 2.*dexc0drs
        END FUNCTION dexc0_ksdt

        DOUBLE PRECISION FUNCTION dexc_ksdt(rs,t,xi)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs,t,xi
          DOUBLE PRECISION :: w,a,b,c,d,e,da,db,dc,dd,de
          DOUBLE PRECISION :: dem,fxc,dfxcdt1,dfxcdt2
          DOUBLE PRECISION :: dfxcdtdrs1,dfxcdtdrs2
          DOUBLE PRECISION :: dfxcdt,dfxcdrs,dfxcdtdrs

          IF (t == 0) THEN
            dexc_ksdt = dexc0_ksdt(rs,xi)
          ELSE
            call coefficients(xi,t,w,a,b,c,d,e,da,db,dc,dd,de)

            dem = 1. + d*sqrt(rs) + e*rs
            fxc = (-1./rs)*(w*a + b*sqrt(rs) + c*rs)/dem
            dfxcdt1 = (-1./rs)*(w*da + db*sqrt(rs) + dc*rs)/dem
            dfxcdt2 = -fxc*(dd*sqrt(rs) + de*rs)/dem
            dfxcdt = dfxcdt1 + dfxcdt2
            dfxcdrs = - (fxc*((d/(2.*sqrt(rs))) + e)/dem)              &
                      - ((1./rs)*fxc)                                  &
                      - ((1./rs)*((b/(2.*sqrt(rs))) + c)/dem)
            dfxcdtdrs1 = - ((1./rs)*(dc + (db/(2.*sqrt(rs))))/dem)     &
                         - (dfxcdt1/rs)                                &
                         - (dfxcdt1*(e + (d/(2.*sqrt(rs))))/dem)
            dfxcdtdrs2 = - (dfxcdrs*(dd*sqrt(rs) + de*rs)/dem)         &
                         - (fxc*((dd/(2.*sqrt(rs))) + de)/dem)         &
                         - (dfxcdt2*((d/(2.*sqrt(rs))) + e)/dem)
            dfxcdtdrs = dfxcdtdrs1 + dfxcdtdrs2

            dexc_ksdt = 2.*(dfxcdrs - t*dfxcdtdrs)
          END IF
        END FUNCTION dexc_ksdt

      END MODULE xc_ksdt_mod
