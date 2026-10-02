!------------------------------------------------------------------------------
! DUMP_SOFT -- validation driver for aao_rad's soft branch (ek < delta).
!
! The radiative branch (sigma(), lines 896 and 1390) is already verified
! pointwise by compare_sigma_stages.py.  The *soft* branch -- lines 776-880 --
! is the other half of the integrand and carries 75% of the cross section in
! the validation run card, so it deserves the same treatment.
!
! The branch is lifted verbatim out of aao_rad.f90 by
! validation/build_dump_soft.sh and pasted into the subroutine below, so there
! is no second copy to drift.  Every input it needs (es, ep, sp, th0, qsq,
! epw, w2, q0, s, csthcm, phicm, ehel, delta) is passed in explicitly, and the
! declaration kinds mirror aao_rad.f90's exactly: es/ep/ps/pp/cst0 are real*8
! (COMMON /radcal/), while q0, s, th0, qsq, epw, w2, sp, csthcm, phicm, f, g,
! fkt, kfac, nu, signr, deltar, delinf and sigr1 are single precision.  That
! matters: sp and w2 are rounded to float32 before they enter the logarithms,
! so a float64 replica would *not* reproduce the original.
!
! The two `go to 20` statements inside the lifted block jump to the reject
! exit that terminates this subroutine, mirroring aao_rad.f90's `go to 20`
! back-edge to the top of the trial loop.  sigr is set to -1 there so the
! caller can tell "rejected" from a genuine zero.
!
! stdin  (whitespace separated, one point per line):
!    es  ep_pre  q2  csthcm  phicm_deg  ehel
!
! stdout (27 columns, ES17.9):
!    1-6   echo of the input
!    7-9   sigr (the branch's own result), sigr1, the value dsigma returned
!   10-15  sigma0 sigu sigt sigl sigi sigip  (maid_lee reassigns these:
!          sigma0=sig_0, sigu=sigma_t, sigt=sigma_tt, sigl=sigma_l,
!          sigi=sigma_lt, sigip=sigma_ltp)
!   16     asym_p from dsigma
!   17-22  f g fkt kfac nu epeps  as the branch computed them
!   23-24  signr, deltar
!   25     delinf
!   26     the log2sp factor, dlog(2 sp / mel^2)
!   27     w2 as the branch sees it (single precision)
!   28     1 if the point was rejected, 0 if it produced a cross section
!
! Must be run from a directory containing spp_tbl/ (or with CLAS_PARMS set),
! because maid_lee() reads the table relative to the working directory.
!------------------------------------------------------------------------------

program dump_soft

   implicit none

   ! Must match aao_rad.f90's declaration order so the COMMON layout agrees.
   COMMON/ALPHA/ ALPHA, PI, MP, MPI, MEL, WG, EPIREA, TH_OPT, RES_OPT

   real*8 alpha, pi, mp, mpi, mel, wg
   integer EPIREA, TH_OPT, RES_OPT

   real*8 es_in, ep_in, q2_in, csthcm_in, phicm_in, ehel_in
   integer ios, n
   real csthcm, phicm
   integer ehel

   real*8 sigr_d, sigr1_d, epeps_d, deltar_d, signr_d, f_d, g_d, fkt_d
   real*8 kfac_d, nu_d, delinf_d, log2sp_d, w2_d
   real sigr, sigr1, epeps, deltar, signr, f, g, fkt, kfac, nu, delinf
   real log2sp, w2, qs
   real sigma0, sigu, sigt, sigl, sigi, sigip, asym_p

   n = 0

   ! --- constants, as aao_rad.f90:232, 290 and mpintp.inc ------------------
   alpha = 1 / 137.d0
   pi = 3.14159d0
   mp = 0.938d0
   mpi = 0.1395d0
   mel = 0.511d-3
   wg = mp + mpi + .0005d0
   EPIREA = 3
   TH_OPT = 7          ! aao_rad.f90:210 reads theory_opt; 7 = MAID07
   RES_OPT = 0         ! aao_rad.f90:218

  10 read(5, *, end=999, iostat=ios) es_in, ep_in, q2_in, csthcm_in&
      , phicm_in, ehel_in

   if (ios.ne.0) goto 999

   n = n + 1
   ehel = nint(ehel_in)
   ! q2 is single precision in aao_rad.f90 (line 84), so narrow it here rather
   ! than inside the wrapper -- the wrapper's `real qs` is already the right
   ! kind, this just matches the two.
   qs = q2_in
   ! csthcm/phicm are single precision in aao_rad.f90 (lines 37 and 77); narrow
   ! them to match the wrapper's `real` dummies.
   csthcm = csthcm_in
   phicm = phicm_in

   call soft_branch(es_in, ep_in, qs, csthcm, phicm, ehel, &
        sigr, sigr1, sigma0, sigu, sigt, sigl, sigi, sigip, asym_p, &
        f, g, fkt, kfac, nu, epeps, signr, deltar, delinf, log2sp, w2)

   sigr_d = sigr
   sigr1_d = sigr1
   epeps_d = epeps
   deltar_d = deltar
   signr_d = signr
   f_d = f
   g_d = g
   fkt_d = fkt
   kfac_d = kfac
   nu_d = nu
   delinf_d = delinf
   log2sp_d = log2sp
   w2_d = w2

   write(6, '(27(1x,es17.9))') es_in, ep_in, q2_in, dble(csthcm_in)&
      , dble(phicm_in), dble(ehel_in) &
      , sigr_d, sigr1_d, dble(sigma0), dble(sigu), dble(sigt), dble(sigl)&
      , dble(sigi), dble(sigip), dble(asym_p) &
      , f_d, g_d, fkt_d, kfac_d, nu_d, epeps_d, signr_d, deltar_d &
      , delinf_d, log2sp_d, w2_d, merge(1.d0, 0.d0, sigr .le. 0.)

   goto 10

 999 continue

   write(0, '(a,i0)') 'dump_soft: read ', n, ' points'

end