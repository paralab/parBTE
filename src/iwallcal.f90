!!$ Ad astra per aspera
!!$ This subroutine calculates the phonon intensity 
!!$ The inputs to this subroutine are polarization, band index and temperature

Subroutine iwalls(pol,mat, band, temp,actint)
!   SUBROUTINE iwalls(ibstep, iband_step, temp,actint)
    Use precisions
    Use constants !,only:zero,two,one,wminla,vsla,four,cla,dirac,half,eight,pi,wminta,vsta,cta,hobol
    Use variables !,only: delta_la,delta_ta,wi,xi,nla, io_debug, debug_level
    Implicit none

    Integer(int_p):: i
    Integer(int_p), Intent(In):: pol,band, mat
    !INTEGER(int_p), INTENT (In):: ibstep, iband_step
  
    !INTEGER(int_p) :: npolar,pol
 
    Real(real_p):: fmin,fmax, a, b, c1, c2, yi, wavenum, disfun, func
    Real(real_p), Intent(In)::temp 
  
    REAL(real_p)  :: delband 
    Real(real_p), INTENT(out):: actint
 
    ! IF(debug_level > 0) WRITE(io_debug,*) "Starting iwallcal_iwall"
    ! nopolar = no_polar (ibstep)
    ! delband = delta_band (ibstep)
    ! actint = zero
    
    ! !!! asssuming only silicon and 2 bandsteps, need to modify for interface
    ! IF (ibstep .eq. 1) THEN
    !     fmin = wminlas + (iband_step-1)*delband
    ! ELSE
    !     fmin = wdisper(ibstep-1) + (band - no_bands_in_step(ibstep-1) -1)*delband  !!! assuming two steps as only silicon
    ! END IF
    fmin = band_low (band)
    fmax = fmin + delta_band(band)
    a = fmin
    b = fmax
    c1 = two/(b-a)
    c2 = one-c1*b
    actint = zero
  
    SELECT CASE (mat)
    CASE (SILICON)
        SELECT CASE (pol)
        CASE (LA_1)
            DO i = 1,20  ! gauss-legendre quadriture
                yi = (xi(i) - c2)/c1 
                wavenum = (( -vslas + sqrt((vslas**2) + four*(yi)*clas))/(two*clas))
                disfun = one/(exp(hobol*yi/(temp)) - one)
                func = yi*disfun*(wavenum**2)   
                actint = actint + wi(i)*func
            ENDDO 
            actint = half*actint*(b-a)
            actint = actint*dirac/(eight*(pi**3))  
        CASE (TA_1)
            DO i = 1,20
                yi = (xi(i) - c2)/c1 
                wavenum = (( -vstas + sqrt((vstas**2) + four*(yi)*ctas))/(two*ctas))
                disfun = one/(exp(hobol*yi/(temp)) - one)
                func = yi*disfun*(wavenum**2)     
                actint = actint + wi(i)*two*func
            ENDDO
            actint = half*actint*(b-a)
            actint = actint*dirac/(eight*(pi**3))
        END SELECT
    END SELECT
End subroutine iwalls     
!  DO i = 1,20  ! gauss-legendre quadriture
!     yi = (xi(i) - c2)/c1
!        DO pol = 1, nopolar
!            IF (pol .ne. LA )        !!! LA and TA only 
!                wavenum = (( -vstas + sqrt((vstas**2) + four*(yi)*ctas))/(two*ctas))
!            ELSE
!                wavenum = (( -vslas + sqrt((vslas**2) + four*(yi)*clas))/(two*clas))
!            ENDIF
!        disfun = one/(exp(hobol*yi/(temp)) - one)
!        func = yi*disfun*(wavenum**2)   
!        actint = actint + wi(i)*func
!        ENDDO
!  ENDDO
!   actint = half*actint*(b-a)
!   actint = actint*dirac/(eight*(pi**3))    
  
!  !Select case (pol)
!  !  case(0) ! la phonons
!        fmin = wminla + (band-1)*delta_la
!        fmax = fmin + delta_la
!        a = fmin
!        b = fmax
!        c1 = two/(b-a)
!        c2 = one-c1*b
!     DO i = 1,20  ! gauss-legendre quadriture
!            yi = (xi(i) - c2)/c1 
!            wavenum = (( -vsla + sqrt((vsla**2) + four*(yi)*cla))/(two*cla))
!            disfun = one/(exp(hobol*yi/(temp)) - one)
!            func = yi*disfun*(wavenum**2)   
!            actint = actint + wi(i)*func
!     ENDDO
!     
!     actint = half*actint*(b-a)
!     actint = actint*dirac/(eight*(pi**3))    

! ! case(1) ! ta phonons
!        fmin = wminta + (band-1-nla)*delta_ta
!        fmax = fmin + delta_ta
!        a = fmin
!        b = fmax
!        c1 = two/(b-a)
!        c2 = one-c1*b     
!     DO i = 1,20
!            yi = (xi(i) - c2)/c1 
!            wavenum = (( -vsta + sqrt((vsta**2) + four*(yi)*cta))/(two*cta))
!            disfun = one/(exp(hobol*yi/(temp)) - one)
!            func = yi*disfun*(wavenum**2)     
!            actint = actint + wi(i)*two*func
!     ENDDO
!     
!     actint = half*actint*(b-a)
!     actint = actint*dirac/(eight*(pi**3))
  
 ! END select
  
 ! ADD other band if neccessary--assy

 !IF(debug_level > 0) WRITE(io_debug,*) "Finishing iwallcal_iwall" 



