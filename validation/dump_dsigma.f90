!------------------------------------------------------------------------------
! DUMP_DSIGMA -- validation driver for the MAID response functions.
!
! The Python port reimplements dsigma()/maid_lee() on the GPU.  This driver
! calls the *original* Fortran dsigma() on a list of kinematic points read from
! stdin and writes the seven response functions plus the intermediate CGLN
! (ff1..ff6) and helicity (hh1..hh6) amplitudes to stdout, so the port can be
! checked point by point instead of only through the sampled event shapes.
!
! stdin  (whitespace separated, one point per line):
!    theta  q2  W  cos(theta_cm)  phi_cm  ehel  opt1  opt2  opt3
!
! stdout (ES16.8 columns):
!    1-5   theta q2 W cos(theta_cm) phi_cm                  (echo of the input)
!    6-12  sig0 sigu sigt sigl sigi sigip asym_p
!          (note: 'sigt' carries sigma_tt and 'sigl' carries sigma_l, because
!           maid_lee.f90 assigns them that way)
!   13-24  Re/Im of ff1..ff6
!   25-36  Re/Im of hh1..hh6
!
! Must be run from a directory containing spp_tbl/ (or with CLAS_PARMS set to
! one), because maid_lee() reads the table relative to the working directory.
!------------------------------------------------------------------------------

program dump_dsigma

   implicit none

   include 'mpintp.inc'
   include 'spp.inc'

   ! Local names cannot collide with the COMMON block variables (q2, W,
   ! csthcm, phicm, epsilon, e_hel, ...), hence the _in suffixes.
   real th_in, q2_in, w_in, cscm_in, phicm_in, ehel_in
   real sig0, sigu, sigt, sigl, sigi, sigip, asym_p
   integer opt1, opt2, opt3, ehel, ios, n

   n = 0
10 read(5, *, end=999, iostat=ios) &
      th_in, q2_in, w_in, cscm_in, phicm_in, ehel_in, opt1, opt2, opt3

   if (ios.ne.0) goto 999

   n = n + 1
   ehel = nint(ehel_in)

   call dsigma(th_in, q2_in, w_in, cscm_in, phicm_in, opt1, opt2, opt3&
      , sig0, sigu, sigt, sigl, sigi, sigip, asym_p, ehel)

   ! The COMMON blocks still hold the state of that last call, so the
   ! intermediates can be recomputed and dumped as well.  That localises a
   ! mismatch to a single stage of the amplitude chain.
   call cgln_amps
   call helicity_amps

   write(6, '(36(1x,es16.8))') th_in, q2_in, w_in, cscm_in, phicm_in&
      , sig0, sigu, sigt, sigl, sigi, sigip, asym_p&
      , dble(real(ff1)), dble(aimag(ff1))&
      , dble(real(ff2)), dble(aimag(ff2))&
      , dble(real(ff3)), dble(aimag(ff3))&
      , dble(real(ff4)), dble(aimag(ff4))&
      , dble(real(ff5)), dble(aimag(ff5))&
      , dble(real(ff6)), dble(aimag(ff6))&
      , dble(real(hh1)), dble(aimag(hh1))&
      , dble(real(hh2)), dble(aimag(hh2))&
      , dble(real(hh3)), dble(aimag(hh3))&
      , dble(real(hh4)), dble(aimag(hh4))&
      , dble(real(hh5)), dble(aimag(hh5))&
      , dble(real(hh6)), dble(aimag(hh6))

   goto 10

999 continue

   write(0, '(a,i0)') 'dump_dsigma: read ', n, ' points'

end
