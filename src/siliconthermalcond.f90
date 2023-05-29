   Subroutine siliconthermalcond(cell,ksi)
 
  Use precisions
  Use constants !only:zero,wminla,two,one,vsla,four,cla,hobol,&
                    ! half,dirac,eight,pi,wminta,vsta,cta
  Use variables !,only:delta_la,delta_ta,xi,wi,nla,io_debug, debug_level
  USE grid
  Implicit None 

  
  INTEGER(int_p), Intent(In):: cell
  INTEGER(int_p):: band, pol, mat
  INTEGER(int_p):: j
  
!  INTEGER(int_p) :: npolar,pol
 
  REAL(real_p):: fmin, fmax, a, b, c1, c2, yi, wavenum, disfun, func, temp, beta
  
  REAL (real_p)  :: delband 
  Real (real_p), intent(out):: ksi
  REAL (real_p) :: bl_l, btn_l, btu_l
  REAL (real_p) :: nor,umk, vel
  REAL (real_p) ::  wmax_half_l, kp_la, kp_ta, k_p_full


  mat = cell_type(cell)
  temp = tnot(cell)

  SELECT CASE (mat) 
  CASE (SILICON)
	  
         bl_l = bl(SILICON)
         btn_l = btn(SILICON)
         btu_l = btu(SILICON)
         wmax_half_l =wmax_half (SILICON)
         pol = LA_1
	    fmin = 0.0
           fmax = wmaxlas 				!fmin + delta_band (band)
           a = fmin
           b = fmax
           c1 = two/(b-a)
           c2 = one-c1*b
	    k_p_full = zero
           Do j = 1,20
                yi = (xi(j) - c2)/c1 
                wavenum = (( -vslas + Sqrt((vslas**2) + four*(yi)*clas))/(two*clas))
                disfun = Exp(hobol*yi/temp)
                vel = sqrt((vslas**2.0) + four * clas *yi)
                beta = bl_l*(temp**3)*(yi**2)  
                func = ((yi*wavenum)**2)*disfun/((disfun-one)**2) *(vel/beta)    
                k_p_full = k_p_full + wi(j)*func
               ! WRITE (*,*) yi, wavenum, disfun, func, k_p
            Enddo
     !STOP
            k_p_full = half*k_p_full*(b-a) 
            k_p_full = k_p_full*(dirac*hobol)/((six*(pi**2))*(temp**2))
	     kp_la = k_p_full
         
  
  
        pol= TA_1
          fmin = 0.0
          fmax = wmaxtas 				!fmin + delta_band (band)
          a = fmin
          b = fmax
          c1 = two/(b-a)
          c2 = one-c1*b
          k_p_full = zero  
            Do j = 1,20
                yi = (xi(j) - c2)/c1 
                wavenum = (( -vstas + Sqrt((vstas**2) + four*(yi)*ctas))/(two*ctas))
                disfun = Exp(hobol*yi/temp)
                vel = sqrt((vstas**2.0) + four * ctas *yi)
                nor   = btn_l*(temp**4)*yi
                umk = zero
                IF (yi .ge. wmax_half_l) Then                  
                    umk = btu_l*(yi**2)/(Sinh(hobol * yi/temp))             
                ENDIF
                beta = nor + umk
                func = ((yi*wavenum)**2)*disfun/((disfun - one)**2)*(vel/beta)
                k_p_full = k_p_full + two*wi(j)*func
            Enddo
            k_p_full = half*k_p_full*(b-a)
            k_p_full = k_p_full*(dirac*hobol)/((six*(pi**2))*(temp**2)) 
   
 
            kp_ta = k_p_full
 	     ksi = kp_la + kp_ta	
     
    CASE (GERMENIUM)
		WRITE (io_debug,'(a)') "COmplete the code"
		STOP
       
     END SELECT  
    
!    IF(debug_level > 0) WRITE(io_debug,*) "Finishing epcalc_ipcal" 

  End Subroutine siliconthermalcond