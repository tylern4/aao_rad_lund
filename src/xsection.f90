subroutine xsection

    implicit none

    include 'mpintp.inc'
    include 'spp.inc'

    real vl, vt, vlt, vtt, vltp
    real phicm_rad, ekin
    real sigma_p, sigma_u

    call multipole_amps ! calc. multipole amplitudes
    call helicity_amps  ! calc. helicity amplitudes

    phicm_rad = phicm * pi / 180.0

    vt = 1.0
    vl = epsilon
    vtt = epsilon
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     vlt and vltp are divided by 2 compared with spp_int_e1
    !     to match AO formalism
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    vlt = sqrt(epsilon * (1 + epsilon) / 2)
    vltp = sqrt(epsilon * (1 - epsilon) / 2)
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ekin = sqrt(Q2) / nu_cm
    fkt = 2 * W * ppi_mag_cm / (W**2 - m_p**2)
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     5 Response fucntions
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    sigma_t = (cabs(hh1)**2 + cabs(hh2)**2&
            + cabs(hh3)**2 + cabs(hh4)**2) / 2.
    sigma_l = cabs(hh5)**2 + cabs(hh6)**2
    sigma_tt = real(hh3 * conjg(hh2) - hh4 * conjg(hh1))
    sigma_lt = sqrt(2.0) * real(conjg(hh5) * (hh1 - hh4) + &
            conjg(hh6) * (hh2 + hh3))
    sigma_ltp = sqrt(2.0) * aimag(conjg(hh5) * (hh4 - hh1) - &
            conjg(hh6) * (hh2 + hh3))
    !
    sigma_l = sigma_l * ekin**2
    sigma_lt = sigma_lt * ekin
    sigma_ltp = sigma_ltp * ekin
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !    sigma_u: unpolarized cross section
    !    simga_h: polarized cross section
    !    asym_ltp: single spin beam asymmetry
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    sigma_u = fkt * (vt * sigma_t&
            + vl * sigma_l&
            + vtt * sigma_tt * cos(2 * phicm_rad)&
            + vlt * sigma_lt * cos(phicm_rad))
    !
    sigma_p = fkt * vltp * sigma_ltp * sin(phicm_rad)
    if (q2.ge.5) sigma_u = sigma_u / (q2 - 4)
    if (q2.ge.5) sigma_p = sigma_p / (q2 - 4)

    sigma_0 = sigma_u + e_hel * sigma_p
    !
    asym_ltp = sigma_p / sigma_u
    !      print *, 'xsection: ',sigma_0,sigma_u, sigma_p,asym_ltp, e_hel
    !
end