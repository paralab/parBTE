SUBROUTINE polarization

    USE precisions,ONLY:real_p,int_p
    USE constants !,only:wmaxlas,wmintas,wmaxtas,wminlas,vslas,vstas,clas,ctas,six,four,half
    USE variables
    
    USE GRID, ONLY: cell_type !...Han added, for checkpolar to identify non-zero system
    !USE aMat, ONLY: m_batchCount
!#ifdef USE_AMAT_BATCHED
    USE amat_batched, ONLY : m_batchCount
!#endif

    IMPLICIT NONE

    INTEGER(int_p) :: band, i , j
    INTEGER(int_p) :: nopolar
    REAL(real_p) :: flow,fhigh,freval
    !  REAL (realp), ALLOCATABLE,DIMENSION(:) :: wdisper
    INTEGER(int_p) ::  no_bands_in_step, pol
    INTEGER(int_p) :: count !...Han added
    
    IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting polarization" 

    ALLOCATE(wdisper(2))          !!! 2 or 4 depends on material if la and ta only considered

    ! wdisper(2)= wmaxtas         !! avoiding the sorting part at the moment.
    ! wdisper(1)= wmaxtag
    wdisper(2)= wmaxlas         !! avoiding the sorting part at the moment.
    wdisper(1)= wmaxtas
    no_b_steps = SIZE(wdisper)
    ALLOCATE(total_bands_in_step(no_b_steps))
    
    !! find band no in each step.
    !total_bands_in_step (:) = zero !...Han comments out b/c it is assigned for each later
    band = 1
    DO i = 1, no_b_steps
        IF (i .eq. 1) THEN
            no_bands_in_step = ((nbands)* (((wdisper(i) - wminlas)) /wdisper(no_b_steps)))  
            total_bands_in_step(i) = no_bands_in_step
        ELSE
            !no_bands_in_step = ((nbands-1)* (((wdisper(i) - wdisper(i-1))) /wdisper(no_b_steps))) +1
            no_bands_in_step = ((nbands)* (((wdisper(i) - wdisper(i-1))) /wdisper(no_b_steps))) +1 !...Siddharth replaced nbands-1 by nband in "memory version"
            total_bands_in_step(i) = no_bands_in_step   !total_bands_in_step (i-1) + no_bands_in_step
        END IF

        !print*,'no_bands_in_step=',no_bands_in_step
    
        !! band discreatization and no of polarization
        DO j = 1, no_bands_in_step
            IF (i .eq. 1) THEN
                delta_band (band) =  (wdisper(i) - wminlas)/ no_bands_in_step
                band_low(band) = wminlas + (band-1)*delta_band(band)
                !print*,'band=',band,',delta_band=',delta_band(band)
            ELSE
                delta_band (band) = (wdisper(i) - wdisper(i-1))/no_bands_in_step
                band_low(band) = band_low(band - 1) + delta_band(band - 1)
            
            ENDIF
            no_polar(band) = no_b_steps + 1 -i
            band = band +1
        ENDDO
    ENDDO
    
    ALLOCATE(polar_type(nbands, np_max))
    
    polar_type (:,:) = zero
    
    DO i = 1, nbands
        nopolar = no_polar(i)
        DO  j = 1, nopolar
            polar_type (i,j) = j
        ENDDO
    END DO
    
    IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing polarization"

    !!!!!!$$$$$!!!! For silicon onlyyy
    !!!$ the input to this subroutine is the number of LA bands NLA and number of TA bands NTA
    !!!!$ assuming that the bins are equally spaced
    !!  delta_la = (wmaxla - wminla)/nla
    !!  delta_ta = (wmaxta - wminta)/nta
    ! ALLOCATE (wdisper(2))          !!! 2 or 4 depends on material if la and ta only considered
    !! ALLOCATE (delta_band(nbands))
    !
    ! wdisper(2) = wmaxlas           !! avoiding the sorting part at the moment.
    ! wdisper(1) = wmaxtas
    ! !WRITE (*,*) "wdisper = ", wdisper
    ! no_b_steps = SIZE(wdisper)
    !! ALLOCATE (delta_band(no_b_steps))
    !! ALLOCATE (no_polar(no_b_steps))
    ! ALLOCATE (total_bands_in_step(no_b_steps))
    ! !WRITE (*,*) "no of steps = ", no_b_steps
    ! !! find band no in each step.\
    ! 
    ! !!! delta_band...no of polarization.
    ! total_bands_in_step (:) = zero
    ! band = 1
    ! DO i = 1, no_b_steps
    !    IF (i .eq. 1) THEN
    !        no_bands_in_step = ((nbands-1)* (((wdisper(i) - wminlas))/wdisper(no_b_steps))) + 1 
    !        total_bands_in_step(i) = no_bands_in_step
    !    ELSE
    !        no_bands_in_step = ((nbands-1)* (((wdisper(i) - wdisper(i-1)))/wdisper(no_b_steps))) + 1
    !        total_bands_in_step(i) = no_bands_in_step   !total_bands_in_step (i-1) + no_bands_in_step
    !        
    !    END IF
    !    !
    !    !! band discreatization and no of polarization
    !    DO j = 1 , no_bands_in_step
    !        IF (i .eq. 1) THEN
    !!            no_bands_in_step(i) = ((nbands-1)* (((wdisper(i) - wminlas))/wdisper(no_b_steps))) + 1 
    !            delta_band (band) =  (wdisper(i) - wminlas)/ no_bands_in_step
    !            band_low(band) = wminlas + (band-1)*delta_band(band)
    !        ELSE
    !!        no_bands_in_step(i) = ((nbands-1)* (((wdisper(i) - wdisper(i-1))/wdisper(no_b_steps))) + 1
    !            delta_band (band) = (wdisper(i) - wdisper(i-1))/no_bands_in_step
    !            band_low(band) = wdisper(i-1) + (band - total_bands_in_step(i-1) -1)*delta_band(band)
    !          
    !        ENDIF
    !            no_polar(band) = no_b_steps + 1 -i
    !            band = band +1
    !    ENDDO
    ! ENDDO
    ! !WRITE (*,*) "total no of band in step" , total_bands_in_step
    ! 
    ! DO i=1,nbands
    !      np_max = MAX(np_max,no_polar(i))
    ! END DO
    ! ALLOCATE (polar_type (nbands, np_max))
    !!WRITE(*,*) "55"
    !!DO i = 1, nbands
    !!    WRITE (*,*) "no polar in bands", i, no_polar(i)
    !!ENDDO
    !!WRITE (*,*) "np_max", np_max
    ! 
    ! polar_type (:,:) = zero
    ! 
    ! !! polar_type(:,1) = LA_1, polar(:, 2) = TA_1   will extend if ge is added.
    !! WRITE (*,*) "bands  nopolar polar_type(i,j)" 
    ! DO i = 1, nbands
    !    nopolar = no_polar(i)
    !    DO  j = 1, nopolar
    !        polar_type (i,j) = j
    !!        WRITE (*,*),  i, j, polar_type(i,j)
    !        
    !    ENDDO
    ! END DO
    !DO i = 1, nbands 
    !     WRITE (*,*),  i,  delta_band(i), band_low(i), no_polar(i)
    !ENDDO
    ! WRITE (*,*),  np_max
    !!  WRITE (*,*) "size of polar type", SIZE(polar_type)

    !...Han added, compute number of systems (to solve for intensity) per band
    !...same as polar_type, all processes compute this for all bands (including bands not owned by me)
    !...this is used in aMat_bached
    !...used throughout the program, so no need to deallocate
    !IF (matMethod == 2) THEN
#ifdef USE_AMAT_BATCHED
    !...number of batches of (non-zero) systems Ax=b per band
    ALLOCATE(nBatchesPerBand(nbands))

    !...values of si and p associated with non-zero systems Ax=b, allocated for max n of non-zero systems
    ALLOCATE(si_nonzeroSystemPerBand(nbands, ndir * np_max))
    ALLOCATE(p_nonzeroSystemPerBand(nbands, ndir * np_max))

    !...determine number of non-zero systems and corresponding (si,p) for each band
    DO band = 1,nbands
        count = 0 !...number of non-zero systems associated with band
        !...loop over all directions, all polarizations to find number of non-zero systems
        DO i = 1,ndir
            DO j = 1,np_max
                pol = polar_type(band, j)
                IF (pol .eq. 0) CYCLE
                !...Oct 11: continue to check the condition of checkpolar
                !...this is only work for the case of 1 material, i.e. if 1 cell has checkpolar(mat,iband,p)=0 then all remaining cells have checkpolar=0
                !...TODO: currently every process check the condition by using cell_type of (local) cell_id=1
                !IF (checkpolar(cell_type(1), band, j) .eq. 0) CYCLE !...comment out on Jan 27, 2022
                count = count + 1 !...increase number of non-zero systems by 1
                si_nonzeroSystemPerBand(band, count) = i !...value of si associated with non-zero system count
                p_nonZeroSystemPerBand(band, count) = j !...value of p associated with non-zero system count
            ENDDO
        ENDDO
        !...check if number of nonzero systems is multiple of "batchSize" (parameter given in variables module)
        IF (mod(count, m_batchCount) .eq. 0) THEN
            !...if yes, then compute the number of batches = (number of nonzero systems / batchSize)
            nBatchesPerBand(band) = count/m_batchCount
        ELSE
            WRITE(*,'(a, i0)') "Error in polarization: # non-zero systems is not divisible by batchSize, band=",band
            STOP
        ENDIF
        ! print*,"band=",band, ", n nonzero systems= ",count, ", n batches=", nBatchesPerBand(band)
        ! if (band .eq. 1) then
        !     do i = 1,count
        !         write(*,'(a,i0,a,i0,a,i0)') "(i,si,p)= ",i,", ",si_nonzeroSystemPerBand(band,i),", ",p_nonZeroSystemPerBand(band,i)
        !     enddo
        ! endif
    ENDDO
!ENDIF
#endif
    ! DO i = 1,nbands
    !     DO j = 1,np_max
    !         write(*,'(a,i0,1x,i0,a,i0,1x,i0)') '(i,j)=',i,j,',checkpolar(:,i,j)=', checkpolar(1,i,j), checkpolar(2,i,j)
    !     enddo
    ! enddo
    !stop
    
    IF(debug_level > 0) WRITE(io_debug,'(a)') "Finishing polarization"           
    !!!!   STOP              
 END SUBROUTINE