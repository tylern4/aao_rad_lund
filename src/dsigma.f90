subroutine dsigma(the, q2, w, cscm, phicm, opt1, opt2, opt3&
   , sig0, sigu, sigt, sigl, sigi, sigip, asym_p, ehel)

   implicit none

   real the, ki_mag, q2, w, cscm, phicm, kf_mag, s2
   real fkt
   real sig0, sigu, sigt, sigl, sigi, sigip, asym_p
   real nu, eps, eps1
   integer ehel
   real cthe
   integer opt1, opt2, opt3
   ! variables for dvmp
   real q2_dummy, W_dummy, hcs_dummy, W2_dummy
   real Mp, mpi0, amp2, amu2, mesonmass, meta, mpip
   real xb_dummy, E_pi_cm, ppi_mag_cm, qv_mag_cm
   real nu_cm, t_dummy
   parameter (Mp = 0.93827)
   PARAMETER (MPI0 = 0.134976, MPIP = 0.13957018, META = 0.547853)

   logical test1, test2, test3
   !
   if (opt2.eq.1) then
      mesonmass = mpi0
   elseif (opt2.eq.3) then
      mesonmass = mpip
   elseif (opt2.eq.5) then
      mesonmass = meta
   endif

   nu = 0.5 * (w**2 + q2 - Mp**2) / Mp
   s2 = sin(0.5 * the)**2
   ki_mag = (nu + sqrt(q2 / s2 + nu**2)) * 0.5

   kf_mag = ki_mag - nu
   eps = 1. / (1 + 2.0 * (1 + nu * nu / q2) * tan(0.5 * the)**2)

   if (kf_mag.lt.0.1) then
      print *, 'low electron momentum', kf_mag, q2, w
      sig0 = 0.
      sigu = 0.
      sigt = 0.
      sigl = 0.
      sigi = 0.
      sigip = 0.
      asym_p = 0.
   endif

   if(opt1.eq.7) then
      call maid_lee(q2, w, eps, cscm, phicm, opt1, opt2, opt3, &
         sig0, sigu, sigt, sigl, sigi, sigip, asym_p, ehel)
   endif
end
