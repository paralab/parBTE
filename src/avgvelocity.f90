!!$ This subroutine calculates the band averaged velocity of each frequency band 

!!! no of polarization in each band is also counted and polar type is store
Subroutine avgvel
    Use precisions,Only:real_p,int_p
    Use constants !,only:wmaxlas,wmintas,wmaxtas,wminlas,vslas,vstas,clas,ctas,six,four,half
    Use variables
    Implicit None

    Integer(int_p)::i, j, band
    INTEGER(int_p) :: pol  !, nopolar, pol
    Real(real_p) :: flow,fhigh,freval
  
    IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting avgvelocity" 

    DO i = 1, nbands
        DO j = 1 , np_max
            pol = polar_type(i,j)
            flow = band_low(i)
            fhigh = flow + delta_band(i)
            freval = (flow + fhigh)*half
            !print*,'(i,j)=',i,j,',flow=',flow,',fhigh=',fhigh,',freval=',freval
            centfreq(i) = freval
            IF (pol .eq. 0) CYCLE
            SELECT CASE (pol)   
                CASE (LA_1)         ! silicon la
                    gpvel(i,LA_1) = sqrt((vslas**2.0) + four * clas * freval)
                CASE (TA_1)         ! silicon ta
                    gpvel(i,TA_1) = sqrt((vstas**2.0) + four * ctas * freval)
            END SELECT
        END DO                     
    END DO         !!! no of bands in step

    wavenumber (:,:) = zero
 
    DO i = 1, nbands
        DO j = 1, np_max
            pol = polar_type(i,j)
            IF (pol == 0) CYCLE
            SELECT CASE (pol)
                CASE (LA_1)
                    wavenumber(i,pol) = (( -vslas + sqrt((vslas**2) + four*centfreq(i)*clas))/(two*clas))
                CASE (TA_1)
                    wavenumber(i,pol) = (( -vstas + sqrt((vstas**2) + four*centfreq(i)*ctas))/(two*ctas))
            END SELECT
        ENDDO
    ENDDO
    !DO i = 1, nbands 
    !     WRITE (*,*),  delta_band(i), band_low(i), centfreq(i)
    !ENDDO
    !STOP
    
    !DO i = 1, nbands
        !    IF (pol .eq. 0) CYCLE
    !    WRITE (*,*), i,rank,  wavenumber(i,LA_1), wavenumber(i,LA_2), wavenumber(i,TA_1),wavenumber(i,TA_2)
    ! ENDDO
    ! 
    ! DO i = 1, nbands
    !    !IF (pol .eq. 0) CYCLE
    !    WRITE (*,*), i,rank,  gpvel(i,LA_1), gpvel(i,TA_1)
    ! ENDDO
    !      
    IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finished avgvelocity" 
    !STOP
End Subroutine avgvel