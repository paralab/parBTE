!***********************************************************************
! Definitions of basic constants

Module constants
  Use precisions,Only:real_p
  Implicit None
  Save

  Real(real_p), Parameter :: zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0D0, & 
                            five=5.0d0, Eight = 8.0d0, nine = 9.0d0, &
                            half=0.5d0,third=0.33333333d0, fourth = 0.25D0, &
                            pi = 3.1415926535898d0, tiny =1.0D-20,tol =1.0d-8, &
                            tzero = 1.0d-15,six = 6.0d0, notsotiny = 1.0d-10

  Real(real_p):: dirac = 1.054571628d-34, boltz = 1.3806503d-23, &
                hobol = 7.63822401661014d-12
                
  ! parameters related to Silicon dispersion from page 39
  Real(real_p),Parameter:: wmintas = 0.0d0, wminlas = 0.0d0, wmaxlas=7.728337675901222d13,& 
                          wmaxtas = 2.97927417405992d13, clas=-2.0d-7,vslas=9.01d3, &
                          ctas = -2.26d-7,vstas= 5.23d3,wmaxlo=9.88d13, &
                          wmaxtos= 10.20d13,wmintos = 8.7124d13, vslos = 0.0d0, &
                          clos = -1.6d-7, vstos = -2.57d3, ctos= 1.11d-7
                          
                                                    
  ! parameters related to Germenium dispersion 
  Real(real_p),Parameter:: wmintag = 0.0d0, wminlag = 0.0d0, wmaxlag=4.41802754820937d13,& 
                          wmaxtag = 1.50052084481175d13, clag=-1.2d-7,vslag=5.30d3, &
                          ctag = -0.82d-7,vstag= 2.26d3,wmaxlog=5.70d13, &
                          wmaxtog= 5.50d13,wmintog = 5.29927d13, vslog = -0.99d3, &
                          clog = -0.48d-7, vstog = -0.18d3, ctog= 0.0d0


                          
  ! relaxation time scale related parameters eqn 2.32, 2.33
  !!wmax_half, grun, lattice_const,rho, btu, bl,btn      !!!!!assuming 2 materials
  Real(real_p), DIMENSION(2):: wmax_half = (/ 2.417d13, 1.005d13/),&
                              lattice_const = (/5.43d-10,5.66d-10/),&
                              rho = (/2.33d3,5.32d3 /),&
                              btu = (/5.5d-18,5.0d-18/), &
                              bl = (/2.0d-24,6.9d-24/),&
                              btn =(/9.3d-13,1.0d-11/)      !!!!!assuming 2 materials
  
  !!  bl(SILICON)= 2.0d-24, btu(SILICON) = 5.5d-18, btn(SILICON) = 9.3d-13, &
  !!  bl(GERMENIUM)= 6.9d-24, btu(GERMENIUM) = 5.0d-18, btn(GERMENIUM) = 1.0d-11,&
  !!  wmax_half(SILICON) = 2.417d13,  lattice_const(SILICON) = 5.43d-10, rho(SILICON) = 2.33d3,&    !!grun(SILICON) = 0.59d0,
  !!  wmax_half(GERMENIUM) = 1.005d13,  lattice_const(GERMENIUM) = 5.66d-10, rho(GERMENIUM) =5.32d3 

  !  Real(real_p), PARAMETER::wmax_half = 2.417d13, bl = 2.0d-24, btu = 5.5d-18, &
  !                           btn = 9.3d-13, grun = 0.59d0, lattice_const = 5.43d-10, &
  !                           rho_s = 2.33d3
                            
  !!! ALUMINUM data (Using Au data)
  ! Real(real_p),Parameter:: G_alum = 2.1d8 ,&   !!interface conductance (3.0d8)
  !                           c_alum = 0.129 , &  !! heat capacity (0.91d3)
  !                           rho_al = 19320      !! density   (2700)
  !...according to Siddharth's code _Implicit_memory
  Real(real_p),Parameter:: G_alum = 2.0d8,&!2.0d8,&!3.0d8,&!2.1d8 ,&   !!interface conductance (3.0d8)
                             c_alum = 115.7d0,&!115.7d0,&!0.91d3,&!0.129 , &  !! heat capacity (0.91d3)
                             rho_al = 18290.8,&!18290.8,&!2700,&!19320      !! density   (2700)
                             k_al = 266.6!266.6!205
  !!! SILICON bulk data  
  Real(real_p),Parameter:: k_si = 148 !! bulk thermal conductivity  

  !!!!!!!!!!!!!!broido data!!!!!!!!!!!!!
  Real(real_p),Parameter:: a_n_las = 7.10d-20, a_u_las = 9.51d-47, a_n_tas = 10.9d-20, a_u_tas = 37.8d-47, debye_temp = 636.0               
    
End Module constants
!***********************************************************************
