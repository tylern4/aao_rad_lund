c     Function defining the true MAID cross section function taking only
c     W (GeV), Q^2 (GeV^2), cos(theta*) and phi* (deg) as arguments.
      
      subroutine maid(beam_energy,w,q2,costheta,phi,helicity,
     *     model_opt,channel_opt,resonance_opt,
     *     sigma0 ,sigma_t, sigma_tt, sigma_l, sigma_lt, sigma_ltp,
     *     asym_p)
      
      implicit none

c     inputs
      real beam_energy, w, q2, costheta, phi

      integer model_opt, channel_opt, resonance_opt, helicity

c     dsigma: ,sig0,sigu,   sigt,     sigl,     sigi,sigip,asym_p)
c     maid: sigma0 ,sigma_t, sigma_tt, sigma_l, sigma_lt, sigma_ltp, asym_p
      
c     outputs
      real sigma0, sigma_t, sigma_l, sigma_tt, sigma_lt,
     +     sigma_ltp, asym_p
      
c     Constants:
      
      real pi
      parameter(pi=4.D0*DATAN(1.D0))
      
      real deg2rad
      parameter(deg2rad=pi/180.0)

c     Variables

      real theta_e, target_mass

c     Functions
      real electron_theta

c     Set target mass
      if(channel_opt.eq.3) then
         target_mass = 939.565e-3
      elseif(channel_opt.eq.2) then
         target_mass = 938.272e-3
      else
         write(*,*) "Error, must specify either channel_opt=2",
     +        "(pi-) or 3 (pi+)"
      endif
      
c     Get electron theta:
      theta_e = electron_theta(beam_energy,target_mass,w,q2)
      
c     Call dsigma:
      call dsigma(theta_e*deg2rad, q2, w, costheta, phi,
     +     model_opt, channel_opt, resonance_opt,
     +     sigma0, sigma_t, sigma_tt, sigma_l, sigma_lt,
     +     sigma_ltp, asym_p,
     +     helicity)
      
      end

      function electron_theta(beam_energy,target_mass,w,q2)
      implicit none

      real electron_theta
      real beam_energy, target_mass, w, q2

c     FUNCTION VARIABLES
      real nu, Eprime, p, pprime, costheta
c     CONSTANTS
      real me
      parameter (me=0.511e-3)
      real pi
      parameter(pi=4.D0*DATAN(1.D0))

      nu = (w**2+q2-target_mass**2)/(2.0*target_mass)

      eprime = beam_energy - nu
c     Debug
c      write(*,*) "E' = ", eprime
c     End debug

      p = sqrt(beam_energy**2 - me**2)
      pprime = sqrt(eprime**2 - me**2)

      costheta = (q2/2 - beam_energy*eprime + me**2)
     +     /(p*pprime)

c     debug
c      write(*,*) "costheta=", costheta
c     end debug

c     Ye's
c      theta = 2*TMath::ASin(TMath::Sqrt((Q2)/(4*Ei*Ef)))*(180/TMath::Pi());
c     Ef=Ei-(w*w+Q2-mn*mn)/(2*mn);
      eprime=beam_energy-(w*w+q2-target_mass*target_mass)
     +     /(2*target_mass);
      electron_theta = 2*asin(sqrt((q2)/(4*beam_energy*eprime)))
     +     * (180/pi)
c     End Ye's
      
c      electron_theta = acos(costheta)
      
      end
