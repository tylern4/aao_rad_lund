!------------------------------------------------------------------------------
! DUMP_SIGMA -- validation driver for aao_rad's sigma() (the Mo & Tsai
! radiative integrand).
!
! compare_sigma_points.py found that the port's motsa_sigma() disagreed with
! the Fortran by a factor that is not a constant, and the only way to localise
! it is to look at the intermediates.  This driver sets up the /ALPHA/ and
! /radcal/ COMMON blocks exactly as aao_rad.f90 does for one trial point, calls
! the *original* sigma(), and then repeats sigma()'s arithmetic here with local
! copies of every intermediate so each one can be compared with the port.
!
! sigma()'s own intermediates are function locals and cannot be dumped, so the
! replica below has to reproduce its return value first; column 12 is that
! check.  Types deliberately mirror sigma()'s declarations (real*4 for ffac,
! gfac, f, g, sig_r, sigf, ...) so a mismatch of order 1e-7 is expected and
! anything larger is a real disagreement.
!
! stdin  (whitespace separated, one trial per line):
!    es  ep_pre  q2  ek  cstk  phik_rad  csthcm  phicm_deg  ehel
! where ep_pre is the scattered-electron energy *before* the exit
! bremsstrahlung (the value sigma() was called with), phik_rad is the photon
! azimuth in radians (the un-shifted value sigma() saw) and phicm_deg is the
! hadronic-frame azimuth in degrees, which is the unit dsigma() takes.
!
! stdout (32 columns, ES17.9):
!    1-9   echo of the input
!   10-11  sigma() (the original, unmodified) and this file's replica of it
!   12-24  qq qsq mf2 epw ffac gfac f g sig_r sigf_pre epeps sigf nu
!   25-31  sig0 sigu sigt sigl sigi sigip asym_p  from dsigma (note the
!          maid_lee.f90 reassignment: sigu=sigma_t, sigt=sigma_tt,
!          sigl=sigma_l, sigi=sigma_lt, sigip=sigma_ltp)
!   32     epsilon that dsigma computed for itself (dsigma.f90:33-34)
!
! Must be run from a directory containing spp_tbl/ (or with CLAS_PARMS set),
! because maid_lee() reads the table relative to the working directory.
!------------------------------------------------------------------------------

program dump_sigma

   implicit none

   ! These must match aao_rad.f90 exactly, in declaration order, so that the
   ! COMMON storage layout agrees.
   COMMON/ALPHA/ ALPHA, PI, MP, MPI, MEL, WG, EPIREA, TH_OPT, RES_OPT
   common /radcal/T0, es, ep, ps, pp, rs, rp, u0, pu, uu, cst0, snt0, csths&
      , csthp, snths, snthp, pdotk, sdotk

   real*8 alpha, pi, mp, mpi, mel, wg, T0
   real*8 es, ep, ps, pp, rs, rp, u0, pu, uu, pdotk, sdotk
   real*8 cst0, snt0, csths, csthp, snths, snthp
   integer EPIREA, TH_OPT, RES_OPT

   ! Kinds mirror aao_rad.f90's declarations so the replica does exactly the
   ! same single/double rounding as the original.
   real*8 es_in, ep_in, q2_in, ek_in, cstk_in, phik_in, Tk
   real cscm_in, phcm_in, ehel_in
   integer ehel, ios, n

   real cthk
   real*8 sp
   real*8 qq, mf2
   real*4 ffac1, ffac2, ffac3, ffac4, ffac5, ffac6, ffac
   real*4 gfac1, gfac2, gfac3, gfac4, gfac
   real*4 f, g, fkt, sig_r
   real epw, q0, kfac, epeps, qsq_d, eps_d, th0_d
   real nu, epeps_r
   real sig0, sigu, sigt, sigl, sigi, sigip, asym_p
   real sigf, sigf_pre, sigr, sigr_re
   real s_kin, s_kin2
   real cstk, sntk, phik

   ! aao_rad.f90's own sigma(), lifted verbatim into the generated file by
   ! validation/build_dump_sigma.sh.  It cannot see the driver's COMMON
   ! declarations, hence the explicit interface.
   interface
      real function sigma(ek, Tk, epcos, epphi, ehel)
      real*8 ek, Tk
      real epcos, epphi
      integer ehel
      end
   end interface

   n = 0

  10 read(5, *, end=999, iostat=ios) es_in, ep_in, q2_in, ek_in, cstk_in&
      , phik_in, cscm_in, phcm_in, ehel_in

   if (ios.ne.0) goto 999

   n = n + 1
   ehel = nint(ehel_in)

   ! --- constants, as aao_rad.f90:232, 290 and mpintp.inc ---------------
   alpha = 1 / 137.d0
   pi = 3.14159d0
   mp = 0.938d0
   mpi = 0.1395d0
   mel = 0.511d-3
   wg = mp + mpi + .0005d0
   EPIREA = 3
   TH_OPT = 7          ! aao_rad.f90:210 reads theory_opt; 7 = MAID07
   RES_OPT = 0         ! aao_rad.f90:218

   ! --- kinematics, as aao_rad.f90:568-656 ------------------------------
   es = es_in
   ep = ep_in
   ps = sqrt(es**2 - mel**2)
   rs = ps / es
   pp = sqrt(ep**2 - mel**2)
   rp = pp / ep
   s_kin = q2_in / 4. / es / ep
   T0 = 2. * asin(sqrt(s_kin))
   snt0 = sin(T0)
   cst0 = cos(T0)
   u0 = es - ep + mp
   pu = sqrt(ps**2 + pp**2 - 2 * ps * pp * cst0)
   uu = u0**2 - pu**2
   csths = (ps - pp * cst0) / pu
   csthp = (ps * cst0 - pp) / pu
   snths = sqrt(1. - csths**2)
   snthp = sqrt(1. - csthp**2)

   cstk = cstk_in
   Tk = acos(cstk)
   sntk = sin(Tk)
   phik = phik_in                    ! already radians, as sigma() expects
   sdotk = es * ek_in - ps * ek_in * cstk * csths&
      - ps * ek_in * sntk * snths * cos(phik)
   pdotk = ep * ek_in - pp * ek_in * cstk * csthp&
      - pp * ek_in * sntk * snthp * cos(phik)

   ! --- the original, untouched ----------------------------------------
   sigr = sigma(ek_in, Tk, cscm_in, phcm_in, ehel)

   ! --- replica, so that every intermediate can be dumped ----------------
   cthk = cos(Tk)
   qq = 2 * mel**2 - 2 * es * ep + 2 * ps * pp * cst0 - 2 * ek_in * (es - ep)&
      + 2 * ek_in * pu * cthk
   mf2 = uu - 2 * ek_in * (u0 - pu * cthk)
   sigr_re = -1.0
   if (mf2 .ge. wg**2 .and. qq .lt. 0.d0) then

      epw = sqrt(mf2)
      sp = es * ep - ps * pp * cst0

      ffac1 = -(mel / pdotk)**2 * (2. * es * (ep + ek_in) + qq / 2)
      ffac2 = -(mel / sdotk)**2 * (2. * ep * (es - ek_in) + qq / 2)
      ffac3 = -2.
      ffac4 = 2 / sdotk / pdotk * (mel**2 * (sp - ek_in**2) + sp * (2 * es * ep - sp + ek_in * (es - ep)))
      ffac5 = (2 * (es * ep + es * ek_in + ep * ep) + qq / 2 - sp - mel**2) / pdotk
      ffac6 = -(2 * (es * ep - ep * ek_in + es * es) + qq / 2 - sp - mel**2) / sdotk
      ffac = ffac1 + ffac2 + ffac3 + ffac4 + ffac5 + ffac6

      gfac1 = mel**2 * (2 * mel**2 + qq) * (1. / (pdotk**2) + 1. / (sdotk**2))
      gfac2 = 4.
      gfac3 = 4. * sp * (sp - 2 * mel**2) / pdotk / sdotk
      gfac4 = (2 * sp + 2 * mel**2 - qq) * (1. / pdotk - 1. / sdotk)
      gfac = gfac1 + gfac2 + gfac3 + gfac4

      qsq_d = -qq
      q0 = es - ep
      th0_d = T0

      call dsigma(th0_d, qsq_d, epw, cscm_in, phcm_in, TH_OPT, EPIREA, RES_OPT&
         , sig0, sigu, sigt, sigl, sigi, sigip, asym_p, ehel)

      s_kin2 = (1. - cst0) / 2.
      epeps = 1 + 2 * (1 + q0**2 / qsq_d) * s_kin2 / (1 - s_kin2)
      if (epeps .le. 1.) then
         epeps_r = -1.
      else
         epeps = 1. / epeps
         epeps_r = epeps
      endif

      kfac = (mf2 - mp**2) / 2. / mp
      nu = (mf2 - mp**2 + qsq_d) / 2. / mp

      f = kfac / (2. * pi**2 * alpha * mp) / (1. + nu**2 / qsq_d) * (sigu + sigl)
      g = mp / (2. * pi**2 * alpha) * kfac * sigu
      fkt = (epw**2 - mp**2 + mpi**2) / 2. / epw
      fkt = sqrt(fkt**2 - mpi**2) * 2. * epw / (epw**2 - mp**2)
      f = f * fkt
      g = g * fkt

      sig_r = ((alpha**3 / (2 * pi * qq)**2) / mp) * (ep / es) * ek_in
      sigf = mp**2 * f * ffac + g * gfac

      sigf_pre = sigf
      if (ffac.gt.0.and.gfac.gt.0.and.sig0.gt.0.and.nu.gt.0.and.sigf.gt.0&
         .and.sigu+epeps*sigl.gt.0.and.epeps.gt.0) then
         sigf = sigf * (1. + (epeps * sigt * cos(phcm_in * pi / 90.)&
            + sqrt(epeps * (1. + epeps) / 2) * sigi * cos(phcm_in * pi / 180.)&
            + ehel * sqrt(epeps * (1. - epeps) / 2) * sigip * sin(phcm_in * pi / 180.))&
            / (sigu + epeps * sigl))
         sigr_re = sig_r * sigf
      endif

      ! epsilon as dsigma computed it for itself (dsigma.f90:33-34); note it
      ! uses Mp = 0.93827, not aao_rad's mp = 0.938.
      eps_d = 1. / (1. + 2.0 * (1. + (0.5 * (epw**2 + qsq_d - 0.93827**2) / 0.93827)**2 / qsq_d) &
         * tan(0.5 * th0_d)**2)
   else
      ffac = -1.; gfac = -1.; f = -1.; g = -1.; sig_r = -1.
      sigf = -1.; sigf_pre = -1.; epeps_r = -1.; nu = -1.; eps_d = -1.
      sig0 = -1.; sigu = -1.; sigt = -1.; sigl = -1.
      sigi = -1.; sigip = -1.; asym_p = -1.
      epw = -1.; qsq_d = -1.
   endif

   write(6, '(32(1x,es17.9))') es_in, ep_in, q2_in, ek_in, cstk_in, phik_in&
      , cscm_in, phcm_in, ehel_in, sigr, sigr_re&
      , dble(qq), dble(qsq_d), dble(mf2), dble(epw) &
      , dble(ffac), dble(gfac), dble(f), dble(g), dble(sig_r) &
      , dble(sigf_pre), dble(epeps_r), dble(sigf), dble(nu) &
      , sig0, sigu, sigt, sigl, sigi, sigip, asym_p, dble(eps_d)

   goto 10

 999 continue

   write(0, '(a,i0)') 'dump_sigma: read ', n, ' points'

end
