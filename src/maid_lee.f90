subroutine maid_lee(a1, a2, a3, a4, a5, theory_opt, channel_opt, &
   resonance_opt, sigma0, sigu, sigt, sigl, sigi, sigip, &
   asym_p, ehel)
   !
   implicit none
   !
   real a1, a2, a3, a4, a5
   !
   include 'mpintp.inc'
   include 'spp.inc'
   !
   integer theory_opt, channel_opt, resonance_opt, n_call
   real E_pi_cm
   real m_pi, alpi
   real sigma0, sigu, sigt, sigl, sigi, sigip, asym_p
   integer ehel
   !
   data n_call /0/
   !
   m_pi = m_pip
   if (channel_opt.eq.1) m_pi = m_pi0
   !
   q2 = a1
   w = a2
   epsilon = a3
   csthcm = a4
   phicm = a5
   !
   if (n_call.ne.0) goto 100
   !
   method_helicity = 1
   !
   write(*,*) 'theory_opt,channel_opt,resonance_opt'
   write(*,*) theory_opt,channel_opt,resonance_opt

   !
   if (theory_opt.eq.7) then
      if (channel_opt.eq.1) then
         call revinm('CLAS_PARMS', &
            'spp_tbl/maid07-PPpi.tbl',data_file)
      endif
      if (channel_opt.eq.2) then
         call revinm('CLAS_PARMS', &
            'spp_tbl/maid07-NPpi.tbl', data_file)
      endif
      if (channel_opt.eq.3) then
         call revinm('CLAS_PARMS', &
            'spp_tbl/maid07-PNpi.tbl',data_file)
      endif
   endif
   !
   !      write(6,*) 'Enter max pi-N angular momentum (0-5): '
   !
   !      read(5,*) mwave_L
   mwave_L = 5
   !
   write(6, *) 'Reading multipoles from ', data_file
   !
   call read_sf_file(data_file, 2)
   !
   n_call = n_call + 1

100 E_pi_cm = 0.5 * (W * W + m_pi**2 - m_p**2) / W
   ppi_mag_cm = E_pi_cm**2 - m_pi**2
   ppi_mag_cm = sqrt(ppi_mag_cm)
   qv_mag_cm = ((W * W + Q2 + m_p**2) / 2.0 / W)**2 - m_p**2
   qv_mag_cm = sqrt(qv_mag_cm)
   !      nu_cm      = sqrt(qv_mag_cm**2-Q2)
   !     this is the  right calculation for nu_cm
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   nu_cm = (W * W - m_p**2 - Q2) / (2 * W)
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   fkt = 2.0 * W * ppi_mag_cm / (W**2 - m_p**2)
   e_hel = ehel
   !
   call xsection
   !
   sigma0 = sigma_0
   sigu = sigma_t
   sigl = sigma_l
   sigt = sigma_tt
   sigi = sigma_lt
   sigip = sigma_ltp
   if (e_hel.gt.-0.5.and.e_hel.lt.0.5) then
      asym_p = 0.0
   else
      asym_p = asym_ltp
   endif
   return
end
