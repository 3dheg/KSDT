KSDT
====

Pade fit to the exchange-correlation energy of the 3D homogeneous electron gas.

## Use:

### Compiling:

Since xc_ksdt_mod.f90 is a Fortran 90 module, compilation is straight-forward:

`gfortran xc_ksdt_mod.f90 YOURCODE.f90 -o YOURCODE`

To compile the example `ksdt.f90`:

`gfortran xc_ksdt_mod.f90 ksdt.f90 -o ksdt`

### Functions:

*DOUBLE PRECISION FUNCTION exc_ksdt(rs,th,xi)*

INPUTS:
* rs = Wigner-Seitz radius, scaled by the Bohr radius
* th = Ratio of temperature to Fermi temperature
* xi = Spin-polarization (1 = polarized, 0 = unpolarized)

RETURNS:
* exc_ksdt = Exchange-correlation energy (Rydbergs)

---

*DOUBLE PRECISION FUNCTION dexc_ksdt(rs,th,xi)*

INPUT:
* rs = Wigner-Seitz radius, scaled by the Bohr radius
* th = Ratio of temperature to Fermi temperature
* xi = Spin-polarization (1 = polarized, 0 = unpolarized)

RETURNS:
* dexc_ksdt = rs derivative of exc_ksdt (Rydbergs/rs)

---

EXAMPLE USAGE:

    PROGRAM main

      USE xc_ksdt_mod
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: rs = 4.0, th = 1.0, xi = 1
      DOUBLE PRECISION :: exc, dexc

      exc = exc_ksdt(rs,t,xi,exc0)
      dexc = dexc_ksdt(rs,t,xi,exc0,dexc0)
      write(*,*) exc, dexc

    END PROGRAM main

## Citation

If you use this module in your calculations, please cite both the fit and original Monte Carlo simulation:

  V. V. Karasiev, T. Sjostrom, J. Dufty, and S. B. Trickey
  [Local Spi-density Approximation Exchange-correlation Free-energy Functional](http://adsabs.harvard.edu/abs/2013arXiv1311.4903K)  
  ArXiv eprints 1311.4903 (physics.chem-ph)  
  Submitted to Phys. Rev. B. Rapid Communications (2013)

  E. W. Brown, B. K. Clark, J. L. DuBois, and D. M. Ceperley  
  [Path Integral Monte Carlo simulation of the warm-dense homogeneous electron gas](http://prl.aps.org/abstract/PRL/v110/i14/e146405)  
  Phys. Rev. Lett. 110, 146405 (2013)

Bibtex:

    @ARTICLE{2013arXiv1311.4903K,
       author = {{Karasiev}, V.~V. and {Sjostrom}, T. and {Dufty}, J. and {Trickey}, S.~B.
            },
        title = "{Local Spi-density Approximation Exchange-correlation Free-energy Functional}",
      journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
       eprint = {1311.4903},
     primaryClass = "physics.chem-ph",
     keywords = {Physics - Chemical Physics, Condensed Matter - Other Condensed Matter},
         year = 2013,
        month = nov,
       adsurl = {http://adsabs.harvard.edu/abs/2013arXiv1311.4903K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @article{PhysRevLett.110.146405,
      title = {Path-Integral Monte~Carlo Simulation of the Warm Dense Homogeneous Electron Gas},
      author = {Brown, Ethan W. and Clark, Bryan K. and DuBois, Jonathan L. and Ceperley, David M.},
      journal = {Phys. Rev. Lett.},
      volume = {110},
      issue = {14},
      pages = {146405},
      numpages = {5},
      year = {2013},
      month = {Apr},
      doi = {10.1103/PhysRevLett.110.146405},
      url = {http://link.aps.org/doi/10.1103/PhysRevLett.110.146405},
      eprint = {1211.6130},
      publisher = {American Physical Society}
    }

