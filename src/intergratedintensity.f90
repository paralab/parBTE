SUBROUTINE intintesity
    USE precisions
    USE constants 
    USE variables
    USE grid
    
    IMPLICIT NONE
    
    INTEGER(int_p) :: idir, iband, ibface,currf, ic_1,ic_2, &
                        p, pol, mat, icell, i, ic_si
    INTEGER(int_p) :: icount
    REAL(real_p) :: jout, jface1, jface2, sdotn1, insdotn1, inot1, temp

    !OPEN(unit=55, file="jface",  status="unknown")
    
    !DO iband = 1, nbands
    DO iband = lowband, highband
        icount = iband - lowband + 1
        DO p = 1, np_max
            pol= polar_type (iband,p)   
            IF (pol == 0) CYCLE
            Do ibface = 1, nbcfaces
                currf = bf_to_f(ibface)
                ic_1 = lfc(currf,1)
                ic_2 = lfc(currf,2)
                
                SELECT CASE (bctype(ibface))   
                CASE(ADIA)
                    jout = zero
                    DO idir = 1,ndir
                        IF(threed) THEN
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf) * ret(idir)
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) + vecfz(currf) * inret(idir)
                        ELSE
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) 
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) 
                        ENDIF      
                        
                        jout = jout + Max(0.0d0,sdotn1)*intensity(ic_1,idir,icount,pol)*(insdotn1/(sdotn1+tiny))
                    ENDDO

                    jfacepos(ibface,pol,iband) = jout
                    jfacepos_s(ibface,pol,iband) = jout

                CASE(SYMM)
                    jout = zero
                    DO idir = 1,ndir
                        IF(threed) THEN
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf) * ret(idir)
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) + vecfz(currf) * inret(idir)
                        ELSE
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) 
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) 
                        ENDIF      
                        
                        jout = jout + Max(0.0d0,sdotn1)*intensity(ic_1,idir,icount,pol)*(insdotn1/(sdotn1+tiny))
                    ENDDO

                    jfacepos(ibface,pol,iband) = jout
                    jfacepos_s(ibface,pol,iband) = jout
        
                CASE(ISOT)
                    jout = zero
                    DO idir = 1,ndir
                        IF(threed) THEN
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf) * ret(idir)
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) + vecfz(currf) * inret(idir)
                        ELSE
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) 
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) 
                        ENDIF                       
                        jout = jout + Max(0.0d0,sdotn1)*intensity(ic_1,idir,icount,pol)*(insdotn1/(sdotn1+tiny))
                    ENDDO

                    jfacepos(ibface,pol,iband) = jout
                    jfacepos_s(ibface,pol,iband) = jout

                CASE (CONJ)
                    jface1 = zero
                    jface2 = zero    
                    DO idir = 1, ndir
                        IF(threed) THEN
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf) * ret(idir)
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) + vecfz(currf) * inret(idir)
                        ELSE
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) 
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) 
                        ENDIF   
                        
                        jface1 = jface1 + MAX(0.0d0,sdotn1)*intensity(ic_1,idir,icount,pol)*(insdotn1/(sdotn1+tiny))
                        jface2 = jface2 - MAX(0.0d0,-sdotn1)*intensity(ic_2,idir,icount,pol)*(insdotn1/(sdotn1+tiny))

                    ENDDO
            
                    jfacepos(ibface,pol,iband) = ABS(jface1)
                    jfaceneg(ibface,pol,iband) = ABS(jface2)
                    jfacepos_s(ibface,pol,iband) = ABS(jface1)
                    jfaceneg_s(ibface,pol,iband) = ABS(jface2)
            
                CASE (HFLUX)
                    jout = zero
                    DO idir = 1,ndir
                        IF(threed) THEN
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf) * ret(idir)
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) + vecfz(currf) * inret(idir)
                        ELSE
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) 
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) 
                        ENDIF   
                            jout = jout + Max(0.0d0,sdotn1)*intensity(ic_1,idir,icount,pol)*(insdotn1/(sdotn1+tiny))
                    ENDDO

                    jfacepos(ibface,pol,iband) = jout
                    jfacepos_s(ibface,pol,iband) = jout

            
                CASE (ALSI)
                    jout = zero
                    ic_si= ic_1
                    DO idir = 1,ndir
                        IF(threed) THEN
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) + vecfz(currf) * ret(idir)
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) + vecfz(currf) * inret(idir)
                        ELSE
                            sdotn1 = vecfx(currf) * rmu(idir) + vecfy(currf) * rxi(idir) 
                            insdotn1 = vecfx(currf) * inrmu(idir) + vecfy(currf) * inrxi(idir) 
                        ENDIF 
                        ! IF (cell_type(ic_1) == ALUMINUM) THEN
                        !    ic_si = ic_2
                        !    sdotn1 = -sdotn1
                        !    insdotn1 = -insdotn1    
                        ! ENDIF  
                        jout = jout + Max(0.0d0,sdotn1)*intensity(ic_si,idir,icount,pol)*(insdotn1/(sdotn1+tiny))
                    ENDDO

                    jfacepos(ibface,pol,iband) = ABS(jout)
                    jfacepos_s(ibface,pol,iband) = ABS(jout)
                    
                END SELECT
            
            END DO 
        END DO
    END DO 
    
    !DO ibface = 1, nbcfaces
    !   DO pol = 1, np_max
    !      DO iband = 1, nbands
    !	   WRITE(55, 201) ibface, pol, iband, jfacepos (ibface,pol,iband), jfaceneg (ibface,pol,iband)	
    !      ENDDO
    !    ENDDO
    ! ENDDO
    !WRITE(55, *) "next"
    ! 201 FORMAT (I6,I6,I6,E15.8,E15.8)
    !!!!!!!! Calculate Inot for everycell, band and pol
    
    DO icell = 1, ncells
        mat = cell_type(icell)
        temp = tnot (icell)
        !DO iband = 1, nbands
        DO iband = lowband, highband
            DO  p = 1, np_max
                pol = polar_type(iband, p) 
                IF (pol == 0)CYCLE
                CALL ipcalc(pol, mat, iband, temp,inot1) 
                inot_c(icell, iband, pol) = inot1
                ! if(icell == 2)then
                !        write(*,*) icell,iband,pol,inot_c(icell,iband,pol)
                !   endif
            END DO
        END DO  
    ENDDO  
    
    ! DO i = 1, nbands
    !    !    IF (pol .eq. 0) CYCLE
    !    WRITE (*,*), i, inot_c(1,i,LA_1), inot_c(1, i,LA_2),inot_c(1,i,TA_1),inot_c(1, i,TA_2)
    !    
    ! ENDDO  
    ! 
    ! DO i = 1, nbands
    !    !    IF (pol .eq. 0) CYCLE
    !    WRITE (*,*), i, inot_c(216,i,LA_1), inot_c(216, i,LA_2),inot_c(216,i,TA_1),inot_c(216, i,TA_2)
    !    
    ! ENDDO  
    !STOP          
END SUBROUTINE intintesity
