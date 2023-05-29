! Subroutine to initialize the gaussian points!  
!***********************************************************************
  Subroutine gauss_points
  Use variables , Only : xi,wi, io_debug, debug_level, rank
  Use precisions,only:int_p
  Implicit None

  Integer(int_p) :: i

 
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting gauss_quad--gauss_points"
  
  ! Uses 20-point Gaussian-Legendre quadrature 

  xi(1) = 0.076526521133497; wi(1) = 0.152753387130725
  xi(2) = 0.227785851141645; wi(2) = 0.149172986472603
  xi(3) = 0.373706088715419; wi(3) = 0.142096109318382
  xi(4) = 0.510867001950827; wi(4) = 0.131688638449176
  xi(5) = 0.636053680726515; wi(5) = 0.118194531961518
  xi(6) = 0.746331906460150; wi(6) = 0.101930119817240
  xi(7) = 0.839116971822218; wi(7) = 0.083276741576704
  xi(8) = 0.912234428251325; wi(8) = 0.062672048334109
  xi(9) = 0.963971927277913; wi(9) = 0.040601429800386
  xi(10) = 0.993128599185094; wi(10) = 0.017614007139152
  Do i = 1,10
     xi(10+i) = -xi(i)
     wi(10+i) = wi(i)
  Enddo
  
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing gauss quad--gauss_points" 

  End Subroutine gauss_points
