Subroutine ipcalc(pol, mat, band, temp,totalint)
 ! SUBROUTINE ipcalc(ibstep, iband_step, temp,totalint)
    Use precisions
    Use constants !only:zero,wminla,two,one,vsla,four,cla,hobol,&
                    ! half,dirac,eight,pi,wminta,vsta,cta
    Use variables !,only:delta_la,delta_ta,xi,wi,nla,io_debug, debug_level
    Implicit None 

    Integer(int_p):: j
    Integer(int_p),intent(in):: band, pol, mat
!   INTEGER(int_p), INTENT (In):: ibstep, iband_step
  
!   INTEGER(int_p) :: npolar,pol
 
    Real(real_p):: fmin, fmax, a, b, c1, c2, yi, wavenum, disfun, func
    Real(real_p), Intent(In)::temp
    REAL (real_p)  :: delband 
    Real(real_p), intent(out):: totalint

!   calculate the total energy  corresponding to that temperature
!   involves integrating over both LA and TA branches

!   integrating over the LA branch
!   IF(debug_level > 0) WRITE(io_debug,*) "Starting epcalc_ipcal" 
    totalint = zero
  
!   nopolar = no_polar (ibstep)
!   delband = delta_band (ibstep)
!   actint = zero
  
!!! asssuming only silicon and 2 bandsteps, need to modify for interface
!   IF (ibstep .eq. 1) THEN
!       fmin = wminlas + (iband_step-1)*delband
!   ELSE
!       fmin = wdisper(ibstep-1) + (band - no_bands_in_step(ibstep-1) -1)*delband  !!! assuming two steps as only silicon
!   END IF
    fmin = band_low (band)
    fmax = fmin + delta_band (band)
    a = fmin
    b = fmax
    c1 = two/(b-a)
    c2 = one-c1*b
    totalint = zero
    SELECT CASE (mat) 
    CASE (SILICON)
        SELECT CASE (pol)
        CASE(LA_1) ! la phonons
            Do j = 1,20
                yi = (xi(j) - c2)/c1 
                wavenum = (( -vslas + Sqrt((vslas**2) + four*(yi)*clas))/(two*clas))
                disfun = Exp(hobol*yi/temp)
                func = ((yi*wavenum)**2)*disfun/((disfun-one)**2)     
                totalint = totalint + wi(j)*func
            Enddo
     
            totalint = half*totalint*(b-a) 
            totalint = totalint*(dirac*hobol)/((eight*(pi**3))*(temp**2))
! integrating over the TA branch - remember TA branch is degenerate  
        CASE(TA_1)   
            Do j = 1,20
                yi = (xi(j) - c2)/c1 
                wavenum = (( -vstas + Sqrt((vstas**2) + four*(yi)*ctas))/(two*ctas))
                disfun = Exp(hobol*yi/temp)
                func = ((yi*wavenum)**2)*disfun/((disfun - one)**2)
                totalint = totalint + two*wi(j)*func
            Enddo
            totalint = half*totalint*(b-a)
            totalint = totalint*(dirac*hobol)/((eight*(pi**3))*(temp**2)) !you might want to store pi^3 value 
        
        END SELECT
    
    END SELECT  
  
End Subroutine ipcalc
