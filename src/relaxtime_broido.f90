!!$ Ad astra per aspera
!!$ this subroutine calculates the band avearges scattering time scales for phonons 
!!$ as per the expressions given by holland
!!$ refer Mazumder and Majumdar JHT 2001, for details

Subroutine relaxtime(pol, mat, band, temp,beta)
  
    Use precisions,Only:int_p,real_p 
    Use constants !,only:wminla,one,bl,wminta,two,btn,zero,wmax_half,hobol,half,btu,third
    Use variables !,only:delta_la,delta_ta,nla,xi,wi,centfreq, io_debug, debug_level
    Implicit None 

    Integer(INT_P),Intent(IN):: pol, band, mat
  
    Real(REAL_P):: nor,umk,freval,mfp_lo,mfp_to
    REAL (REAL_P) :: bl_l, btn_l, btu_l, wmax_half_l
    Real(REAL_P),Intent(IN):: temp
    Real(REAL_P),Intent(OUT):: beta 
  
!   IF(debug_level > 0) WRITE(io_debug,*) "Starting relaxtime" 
 
    beta = -1

    SELECT CASE (mat)
        CASE (SILICON)
            bl_l = bl(SILICON)
            btn_l = btn(SILICON)
            btu_l = btu(SILICON)
            wmax_half_l =wmax_half (SILICON)
            SELECT CASE(pol) 
            
            CASE(LA_1) !!$ LA phonon SI
                freval = centfreq(band)      
!               beta = bl_l*(temp**3)*(freval**2) 
                beta = a_n_las*freval**2*temp*(1 - exp(-3*temp/debye_temp)) &
                                + a_u_las*freval**4*temp*(1 - exp(-3*temp/debye_temp)) 
            CASE(TA_1)  !!$ TA phononSi 
                freval = centfreq(band) 
!               nor   = btn_l*(temp**4)*freval
!               umk = zero
!               IF (freval .ge. wmax_half_l) Then                  
!                   umk = btu_l*(freval**2)/(Sinh(hobol * freval/temp))             
!               ENDIF
!               beta = nor + umk
                beta = a_n_tas*freval**2*temp*(1 - exp(-3*temp/debye_temp)) &
                                + a_u_tas*freval**4*temp*(1 - exp(-3*temp/debye_temp)) 
                
            END SELECT
           !beta=1.0
    END SELECT
        


End Subroutine relaxtime
