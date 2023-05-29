!!$**********************************************************************
!!$**********************************************************************
!!$           Main Program to initiate control angle discrete ordinate 
!!$           quadrature scheme
!!$           Fluids and Thermal Analysis Laboratory
!!$**********************************************************************
Subroutine dircalc
  
  Use precisions
  Use variables      
  Use constants,Only: half,two,pi,zero,tzero,one
  USE grid
  Implicit None

  Integer(int_p) :: i, count, j
  
  Real(real_p) :: theta,anphi
  
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting snweights_dircalc" 

  !  ndir = Int(Ntheta*Nphi) ! variable declaration??
  
  !...OPEN(unit=55, file="unit_sand_nb.dat",  status="unknown")
  
  !  Allocate (Omega(ndir))
  !  Allocate(rmu(ndir),rxi(ndir),ret(ndir))
  !  Allocate(inrmu(ndir),inrxi(ndir),inret(ndir))    
  !  dphi = two*pi/Nphi
  !  dtheta = pi/Ntheta 

  count = 0
  
  Do i = 1, Ntheta
    theta = (i-half)*dtheta     
    Do j = 1,Nphi
      anphi = (j-half)*dphi
      count = count + 1
      Omega(count) = two*Sin(theta)*Sin(half*dtheta)*dphi        
      rmu(count) = Sin(theta)*Sin(anphi)
      rxi(count) = Sin(theta)*Cos(anphi)
      ret(count) = Cos(theta)
      inrmu(count) = Sin(anphi) * Sin(dphi*half) * (dtheta - Cos(two*theta)*Sin(dtheta))
      inrxi(count) = Cos(anphi) * Sin(dphi*half) * (dtheta - Cos(two*theta)*Sin(dtheta))
      inret(count) = half * dphi * Sin(two*theta) * Sin(dtheta)
    Enddo
  Enddo

  !...Han added: precompute sdotn and insdotn
  ! DO i = 1,nfaces
  !   DO j = 1,ndir
  !     IF (threed) THEN
  !       sdotn_array(i,j) = (vecfx(i) * rmu(j)) + (vecfy(i) * rxi(j)) + (vecfz(i) * ret(j))
  !       insdotn_array(i,j) = (vecfx(i) * inrmu(j)) + (vecfy(i) * inrxi(j)) + (vecfz(i) * inret(j))
  !     ELSE 
  !       sdotn_array(i,j) = (vecfx(i) * rmu(j)) + (vecfy(i) * rxi(j))     
  !       insdotn_array(i,j) = (vecfx(i) * inrmu(j)) + (vecfy(i) * inrxi(j))
  !     ENDIF
  !   ENDDO
  ! ENDDO
  
  !  DO i = 1, ncells
  !    WRITE (55,*) xc(i), yc(i), zc(i)
  !!  ENDDO
  !  WRITE (55,*) "vector n, nfaces", nfaces
  !  DO i=1,nfaces
  !        WRITE(55,*) vecfx(i),vecfy(i),vecfz(i)
  !  END DO
  !   WRITE (55,*) "s_vector, dir", ndir

  ! DO i =1 , ndir
  !   WRITE(55,*) i, rmu(i),rxi(i),ret(i), omega(i)
  ! ENDDO
    
  !  CLOSE(55)
  
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing snweights_dircal"
  
  ! STOP

End Subroutine dircalc