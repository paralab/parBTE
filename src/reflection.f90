  Subroutine reflec(sin,ibface,sdotn,sout)

  Use precisions
  Use constants,only:two, one
  Use grid !,only:vecfx,vecfy,vecfz
  Use variables !,only:rmu,rxi,ndir,rxi,ret, io_debug, debug_level,one 

  Implicit None 

  INTEGER(int_p)::i,sout
  INTEGER(int_p), INTENT(in) :: sin,ibface
  
  REAL(real_p) :: srx,sry,srz, error, sdots0
  REAL(real_p), INTENT(in):: sdotn
  
  
  !IF(debug_level > 0) WRITE(io_debug,*) "Starting reflection_reflec" 

  ! now the components of intensity are known the value of G( in this case theta_b) at cell centers needs to be calculated   ! calculation of G
  srx = rmu(sin) - two*sdotn*vecfx(ibface)
  sry = rxi(sin) - two*sdotn*vecfy(ibface)
 ! srz = ret(sin) - two*sdotn*vecfz(ibface)
  IF (threed) THEN 
  	srz = ret(sin) - two*sdotn*vecfz(ibface)
  ENDIF
  sout = 0
  error = 1e6
  do i = 1, ndir
  !    IF (threed) THEN
  !       IF ((ABS(srx - rmu(i)) .lt. 1.0d-6) .and. (ABS(sry - rxi(i)).lt. 1.0d-6) .and. (ABS(srz - ret(i)).lt. 1.0d-6) ) then
  !              sout = i
  !              EXIT       
  !       ENDIF
  !  ELSE
  !       IF ((ABS(srx - rmu(i)) .lt. 1.0d-6) .and. (ABS(sry - rxi(i)).lt. 1.0d-6) ) THEN
  !              sout = i
  !              EXIT       
  !       ENDIF
  !  ENDIF
     IF (threed) THEN
     	sdots0 = rmu(i)*srx + rxi(i)*sry + ret(i)*srz
     ELSE
       sdots0 = rmu(i)*srx + rxi(i)*sry
     ENDIF
     IF (ABS(sdots0 - one) .lt. error) THEN
        sout = i
        error = ABS(sdots0 - one)
   ENDIF
!!! fixing reflection.f90

  end do
  !WRITE (*,*) "ibface, si, sout", ibface, sin, sout
  if (sout==0) then
        WRITE (io_debug,'(a)') "sout = 0"
        stop
  endif
  
 ! IF(debug_level > 0) WRITE(io_debug,*) "Finishing reflection_reflec" 
 
  END subroutine reflec
