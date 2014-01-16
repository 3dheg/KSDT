      PROGRAM main

         USE xc_ksdt_mod
         IMPLICIT NONE
         DOUBLE PRECISION, PARAMETER :: rs = 4.0, t = 0.0, xi = 0
         DOUBLE PRECISION :: exc, dexc

         exc = exc_ksdt(rs,t,xi)
         dexc = dexc_ksdt(rs,t,xi)
         write(*,*) exc, dexc

      END PROGRAM main
