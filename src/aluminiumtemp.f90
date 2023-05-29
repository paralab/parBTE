!! subroutine to calculate boundary temp for vitual aluminum to simulate laser heating!
!!$ The inputs to this subroutine are boundary face

    Subroutine aluminiumtemp(ibface)
!   SUBROUTINE iwalls(ibstep, iband_step, temp,actint)
    Use precisions
    Use constants !,only:zero,two,one,wminla,vsla,four,cla,dirac,half,eight,pi,wminta,vsta,cta,hobol
    Use variables !,only: delta_la,delta_ta,wi,xi,nla, io_debug, debug_level
    Use grid
  
    Implicit none

    Integer(int_p):: i,j
    Integer(int_p), Intent(In):: ibface
    INTEGER(int_p):: iband, p  
    INTEGER(int_p) :: currf, ic, ic_1, ic_2, curr_f
    Real(real_p):: t, V, c, A_flux, a_face, dely,  resis, a ,b, temp, G
    REAL (real_p)  :: t_i, t_c, radius, t_sub 
    Real(real_p):: k_silicon !, qnot
    REAL(real_p) :: fluxnet
  
    REAL(real_p) :: tempe
 
    IF(debug_level > 0) WRITE(io_debug,*) "Starting alunimum temp"
  
    currf = bf_to_f(ibface)
    ic_1 = lfc(currf,1)
 

!   ic_2 = lfc(currf, 2)
!   IF (cell_type(ic_1)== SILICON) THEN
        ic_1 = lfc(currf,2)
!       ic_2 = lfc(currf,1)
!   ENDIF
  
  
    t = dt
!   rho_al = 2700
!   10e-9 * 30e-6 * 1     !!! needs to be taken from model file
    a_face= areaf(currf)  !30e-6; !% area_out
    V = a_face * layer 		!!!! 100d-9  ! layer !80d-9  ! changed from 10d-9  !!!! volume of aluminum with 100nm hieght 
!   WRITE (*,*) "Vol", V
!   L = 300e-6;
!   c= 0.91e3 !!! heat capacity of aluminium, add it to constant module
!   DO j = 1,nfcell(ic_1)
!       curr_f =lcf(ic_1,j)
!       IF (xf(curr_f) .lt. 1.0d-9) THEN
!           A_flux = areaf(curr_f)
!           EXIT 
!       ENDIF     
!   ENDDO
!   A_flux =  30e-6 * 1           !! from model file, value needs attention when doing actual tdtr
    qnot = qdblprime
    radius = yf(currf)
   
!   IF(threed) radius = SQRT (yf(currf)**2 + zf(currf)**2)
    IF(threed) radius = SQRT (xf(currf)**2 + yf(currf)**2 + zf(currf)**2)  !!! assuming that the alsi plane remains on a zero plane

  
    IF (radius .gt. spotsize) qnot = zero 
    !WRITE (*,*) ibface, radius, qnot
 
!   WRITE (*,*) "a_face", a_face 
!!!!!!!!!% resistance !!!!!

!   G = 3e8 !% conductance, 110  MW/m2K , add to constant module 
!   dely =10e-9 !; %25e-6;
!   k_si = 148 !!!%W/mK, bulk value, constant module or calculated value???
!   WRITE (*,*) ibface, "alum1"
!   currf = bf_to_f(ibface)
!   ic = lfc(currf,1)
!   IF (cell_type(ic)== ALUMINUM) ic = lfc(currf,2)
  
!   WRITE (*,*) ic, currf, "alum2"
  
    dely = ABS(xf(currf)- xc(ic_1))     
!   WRITE (*,*) "dely", dely
 
!!  calcualtion of silicon thermal condtuctivity on fly
!   k_silicon = k_sil(ic_1)
    k_silicon = k_si

!   resis = (1/(G_alum*a_face)) + ((dely)/(k_silicon*a_face))
    resis = (1/(G_alum*a_face))
!   WRITE (*,*) ibface, ic_1 , k_silicon, resis, a_face, rank

    a= 1/(rho_al*V*c_alum*resis);
    qnot=qnot*exp(-2.0*radius*radius/(spotsize*spotsize)) !...Siddharth added on Aug 8, 2022, for Gaussian Laser Flux
    b= (qnot*a_face)/(rho_al*V*c_alum);

    t_c = tnot(ic_1);
    t_sub = temp_bc(ibface)
!   t_i = temp_bcone(ibface);       !! previous timestep value
    t_i = ttopone(ibface)
!   WRITE(*,*) t_c, t_i

    temp   = t_sub + ((t_i-t_sub)*(exp (-a*t)) + (b/a)* (1 - (exp (-a*t)))); 
!   temp   = t_c + ((t_i-t_c)*(exp (-a*t)) + (b/a)* (1 - (exp (-a*t))));
 
!   temp_bc(ibface) = temp
    ttop(ibface) = temp       !!! ttop is the tranducer!!!!
 
!   fluxnet = (ttop(ibface) - t_c)/(resis*a_face)
!!! flux from jout and jin

    fluxnet = zero
    Do iband = 1, nbands
        DO p = 1, np_max
            fluxnet = fluxnet - flux(ibface,p,iband)
            ! if (ibface .eq. 11) then
            !     write(*,*) '(iface,p,iband,flux())= ', ibface,p,iband, flux(ibface,p,iband)
            !     write(*,*) 'fluxnet= ',fluxnet
            ! endif
        ENDDO
    ENDDO
!   IF (ibface .eq. 10) PRINT *, "fluxnet", fluxnet
!   fluxnet = ABS(fluxnet)
!   fluxnet = - qnot
!!! ttop is the transducer temp
!!! temp_bc is the substrate boundary temp
!!! tbot == temp_bc

    !WRITE(*,*) ibface, temp, t_c, temp_bc(ibface)	
    tempe = ttop(ibface) -(fluxnet/(G_alum))
    temp_bc(ibface) = tempe *rel +  temp_bc(ibface)*(1-rel)
    !write(*,*) '(ibface, fluxnet)= ',ibface, fluxnet, temp_bc(ibface)
!   Using temp as input
!   In the statement below qdblprime is temp input
!   temp_bc(ibface)=qdblprime 

!   ttop(ibface) = temp_bc(ibface) - ((qnot)/(G_alum + (rho_al * c_alum* 80e-9*(1/t))))
      

 
   
!   WRITE(*,*) fluxnet, ibface, dely, k_silicon, a_face, resis
!    WRITE(*,*) fluxnet,ibface,qnot,ttop(ibface),temp_bc(ibface),tempe
!   tbot(ibface) = t_c + fluxnet *(dely/k_silicon)

  
    IF(debug_level > 0) WRITE(io_debug,*) "Finishing alunimum temp"

  
    End subroutine aluminiumtemp
