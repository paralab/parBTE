SUBROUTINE tempcalc
  
  Use precisions
  Use constants,Only:zero,four,tzero,pi,two,one
  Use grid
  Use variables

  Implicit None   
  
  INTEGER(int_p) :: ic, iband, pol, j, mat, f, i, ic_1, ic_2

  INTEGER (int_p) :: ibface, bcon, currf, ifc
  INTEGER (int_p) :: itera, niter
  
  REAL(real_p):: guess, inot, iprime, corr, imfp, ilp
  REAL(real_p):: vel, beta, gnb, gna
  REAL(real_p):: told, uold, unew, uprime, inotone
  
  !!! for temp_bc
  REAL(real_p) :: jin, ftemp, qnot_c
  REAL(real_p) :: actin, actinprime
 
  REAL(real_p):: tempmax, tempmin	

  !CHARACTER(len=6) :: nProcsName, rankName
  !WRITE(nProcsName, '(I6)') numprocs
  !WRITE(rankName, '(I6)') rank

  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting tempcalc" 
  
  !...OPEN(unit=88, file="newton.dat",  status="unknown")

  !OPEN(unit=89, file="t.dat",  status="unknown")
  !OPEN(unit=90, file="gc_band.dat",  status="unknown")
  !OPEN(unit=91, file="in_band.dat",  status="unknown")

  !...update number of iterations in Newton (nonlinear) solving for temperature
  iter_cal = iter_cal + 1

  !...WRITE(88,*) "rank= ", rank ,", time = ", timeindex, ", iter = ", iter_cal

  !...preset max number of iterations
  niter = 100

  !open(32,file='tempcalc_'//TRIM(adjustl(nProcsName))//'_'//TRIM(adjustl(rankName))//'.txt',status = 'unknown')

  !...loop over cells owned by my rank, accumulate residual from each cell
  residual = zero
  DO ic = 1,ncells
    !guess = tnot(ic)
    mat = cell_type(ic) 
    ! IF  (mat== ALUMINUM) CYCLE         !! material_type
    told  = tnotone(ic)     
    
    !gna = zero     
    gnb = zero     
    uold = zero  
    !unew = zero
    !uprime = zero

    DO iband = 1,nbands
      DO j = 1,np_max 
        pol = polar_type(iband, j)
        IF (pol .eq. 0) CYCLE
        !f =  checkpolar(mat, iband, j)
        !IF (f /=0 )THEN
          vel = gpvel(iband, pol)
          !CALL relaxtime(pol, mat, iband, guess, beta)
          ! beta = 1.0d-6
          !CALL iwalls(pol, mat,iband, guess,inot)  
          CALL iwalls(pol, mat,  iband, told, inotone)  
          !CALL ipcalc(pol,mat, iband, guess,iprime)     
          !imfp = beta/vel
          ilp  = idt/vel
          !gna = gna + gc(ic,pol,iband)*(imfp + ilp)
          gnb = gnb + gcone(ic, pol, iband) * ilp
          uold = uold + (four * pi * ilp * inotone)
          !unew = unew + four*pi*(imfp+ilp)*inot
          !uprime = uprime + four*pi*(imfp+ilp)*iprime
        !ENDIF    
      END DO ! polarization loop    
    ENDDO ! band loop
    !write(*,'(a,2(e20.10,1x))') 'printing in tempcalc: (uold,gnb)=',uold,gnb
    !!!!!!!!!!!!!!!!! newton iteration
          
    DO itera = 1,niter
      guess = tnot(ic)
      gna = zero
      unew = zero
      uprime = zero
      DO iband = 1, nbands
        DO j = 1, np_max 
          pol = polar_type(iband, j)
          IF (pol .eq. 0) CYCLE
          !f =  checkpolar (mat, iband, j)
          !IF (f /=0 )THEN
            vel = gpvel(iband, pol)
            CALL relaxtime(pol, mat, iband, guess, beta)
            imfp = beta/vel
            ilp  = idt/vel  
            !CALL iwalls(pol, mat,iband, guess,inot)
            CALL ipcalc(pol, mat, iband, guess, iprime)    
            CALL iwalls(pol, mat, iband, guess, inot)
            !write(*,'(a,2(e20.10,1x))') 'printing in tempcalc: (gna,gc)=',gna,gc(ic, pol, iband)
            gna = gna + gc(ic, pol, iband)*(imfp + ilp)
            ! if (ic==1) then
            !   write(*,'(a,3(i3,1x),2(e20.12))') '(itera,iband,pol,gc,gna)=',itera,iband,pol,gc(ic,pol,iband),gna
            ! endif
            unew = unew + four * pi * (imfp + ilp) * inot
            uprime = uprime + four * pi * (imfp + ilp) * iprime
          !ENDIF
        ENDDO
      ENDDO
      !write(*,'(a,2(e20.10,1x))') 'printing in tempcalc: (unew,gna)=',unew,gna
      !... this is correction of temperature
      corr = (uold + gna - gnb - unew)/uprime
      
      !corr = (uold/uprime) + (gna/uprime) - (gnb/uprime) - (unew/uprime)
      tprime(ic) = corr
      tnot(ic) = tnot(ic) + (relax_temp * tprime(ic))

      ! IF (lc_to_gc(ic) == 1) THEN 
      !   WRITE(*,'(a,/,3(i3,1x),5(e20.12,1x))'), '(rank, global_cell_id, itera, uold, gna, gnb, unew, uprime):', &
      !     rank, lc_to_gc(ic), itera, uold, gna, gnb, unew, uprime
      !   !WRITE(88,*) "cell temp,", itera,",",tprime(ic),",", tnot(ic)
      ! ENDIF

      !...condition to exit itera loop
      IF (ABS (corr) .lt. 1.0d-4) THEN
        !if (rank .eq. 0) print*,'cell, itera=',ic,itera
        EXIT
      ENDIF

    ENDDO !...itera = 1,nitera
    
    !...accumulate residual for cell ic (this residual is used for checking convergence of outer loop)
    residual  = residual + corr**2

    !write(32,'(4(i3,1x),2(e15.6,1x))') rank,ic,lc_to_gc(ic),itera,corr,residual

    !tprime(ic) = corr  
    !tnot(ic) = tnot(ic) + relax_temp*tprime(ic)    

  ENDDO !... ic=1,ncells

  !!!! temp_bc

	tbprime(:) = zero

	DO ibface = 1,nbcfaces
    bcon = bctype(ibface)
    currf = bf_to_f(ibface)
    ic = lfc(currf,1)
    !ftemp = temp_bc(ibface)
    mat = cell_type (ic)
    !PRINT *, "in temp bc cal" ,ibface, bcon
    !IF (bcon == ADIA .and. mat== ALUMINUM) CYCLE 
    IF ((bcon == ISOT) .or. (bcon == CONJ) .or. (bcon == SYMM) ) CYCLE
    IF (bcon ==ALSI)  THEN 
      CALL aluminiumtemp (ibface) !...right now, this function is never called due to BCs is not ALSI (alumnium silicon)
      ! print *, "aluninum temp pass"		
      CYCLE
    ENDIF
    jin = zero
    !actin = zero
    !actinprime = zero
    
    !IF(bcon == ADIA) qnot = zero  
    qnot_c = zero    
    DO iband = 1, nbands
      DO j = 1, np_max 
        pol = polar_type(iband, j)
        IF (pol .eq. 0) CYCLE
        jin = jin +  jfacepos (ibface,pol,iband) 
        !CALL iwalls(pol,mat, iband,ftemp,inot)
        !actin = actin + inot
        !CALL ipcalc(pol, mat, iband, ftemp,iprime)
        !actinprime = actinprime + iprime
      ENDDO        !! polarization loop
    ENDDO       !!! band loop
    DO itera = 1, niter
      ftemp = temp_bc(ibface)
      actin = zero
      actinprime = zero
      DO iband = 1, nbands
        DO j = 1, np_max 
          pol = polar_type(iband, j)
          IF (pol .eq. 0) CYCLE
          CALL iwalls(pol,mat, iband,ftemp,inot)
          actin = actin + inot
          CALL ipcalc(pol, mat, iband, ftemp,iprime)
          actinprime = actinprime + iprime
        ENDDO        !! polarization loop
      ENDDO
      !corr = (qnot  + jin - (pi* actin))/ (pi* actinprime)
      corr = (jin + qnot_c  - (pi* actin))/ (pi* actinprime)
      tbprime(ibface) = corr
      temp_bc(ibface) = temp_bc(ibface) + relax_temp*tbprime(ibface)
      ! IF (ibface == 10) THEN
      !   !WRITE(*,*) "ibface",ibface, itera,tbprime(ibface), temp_bc(ic)
      !   WRITE(88,*) "tbc,", itera,",",tbprime(ibface),",", temp_bc(ibface)
      ! ENDIF
      IF (ABS (corr) .lt. 1.0d-4) EXIT
    ENDDO
    !residual  = residual + corr**2   !...Han? who commented this line?
 	ENDDO !! bface loop
  
  !...Aug 11, 2021: bug found here, taking square root of residual before reduction make the total residual
  !...be bigger than sequential case (i.e. sqrt(r) > r if r < 1, then sqrt(r1) + sqrt(2) > sqrt(r1 + r2)
  !...residual = SQRT(residual)

  !IF (rank.eq.0) print*,'iter_call=',iter_cal,', residual=',residual
  !print"('tempcalc:rank=',i3,', residual sum over cells=',e20.8)",rank,residual

  tempmin = 50
  tempmax = 600
  
  DO i = 1,nbcfaces
    IF(bctype(i) .eq. SYMM) THEN
      ic = bf_to_c(i)
      ifc = bf_to_f(i)
      temp_bc(i) = tnot(ic)
      temp_bc(i) = MAX(tempmin,MIN(tempmax,temp_bc(i)))
    ENDIF
    IF (bctype(i) .eq. CONJ) THEN
      ifc = bf_to_f(i)
      ic_1 = lfc(ifc,1)
      ic_2 = lfc(ifc,2)
      temp_bc(i) = wcf(ifc)* tnot(ic_2) + (one-wcf(ifc))*tnot(ic_1)
    ENDIF
  ENDDO
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing tempcalc" 
  202 FORMAT (I6, I6, E15.8,E15.8,E15.8,E15.8)
  203 FORMAT (I6, I6,I6, E15.8,E15.8,E15.8,E15.8)

END SUBROUTINE tempcalc