      MODULE constants
          ! Global Constants
          integer, parameter :: MAX_INTEGER_DIGITS = 10
      END MODULE constants

      PROGRAM main

         USE xc_ksdt_mod
         USE constants           ! Constants defined above

         ! Disable implicit declarations (i-n rule)
         IMPLICIT NONE
         DOUBLE PRECISION :: rs, t, xi
         DOUBLE PRECISION :: exc, dexc

         ! Variable defintions
         CHARACTER(MAX_INTEGER_DIGITS) :: rs_str, t_str, xi_str

         ! First command line argument is the base, second is the exponent
         call getarg(1, rs_str)
         call getarg(2, t_str)
         call getarg(3, xi_str)

         ! Make sure user provided both base and exponent
         if (rs_str=='' .OR. t_str=='' .OR. xi_str=='') then
             stop 'Usage: rs t xi'
         endif

         ! Convert strings to integers
         read (rs_str, *) rs
         read (t_str, *) t
         read (xi_str, *) xi

         exc = exc_ksdt(rs,t,xi)
         dexc = dexc_ksdt(rs,t,xi)
         write(*,*) exc, dexc

      END PROGRAM main
