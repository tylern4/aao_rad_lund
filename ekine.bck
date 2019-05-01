      program test
c
      implicit none
c
c
      real ki_mag, kf_mag, theta_e
      real m_p,dummy1, dummy2, W, Q2
      integer i
c
      m_p  = 0.93827
c
      write(6,*) 'input electron angle'
      read(5,*) theta_e
      ki_mag = 6.045
      theta_e = 3.141592*theta_e/180.
c
      do i = 1, 15
         W = 1.6 + (i-1) * 0.1
         dummy1 = 2*m_p*ki_mag - m_p**2 - W**2
         dummy2 = 2*m_p*(1+1/(4*ki_mag*sin(theta_e/2)**2))
         Q2 = dummy1 / dummy2
         kf_mag = Q2 / (4*ki_mag*sin(theta_e/2)**2)
         print *, 'W= ',W,' Q2= ',Q2 , ' kf_mag= ',kf_mag
      enddo
c
      stop
      end
