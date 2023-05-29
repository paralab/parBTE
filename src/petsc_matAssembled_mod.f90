! \brief module to solve using matrix-assembled PETSc method
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! Assemble PETSc matrix and RHS vector to solve the linear systems of intensity
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
MODULE petsc_matAssembled
    USE petsc_solver
    IMPLICIT NONE

    CONTAINS

    !...allocate, set type and size of the matrix
    SUBROUTINE allocate_petsc_mat
        IMPLICIT NONE

        CALL determine_max_faces_per_cell(m_maxFacesPerCell)

        CALL MatCreate(m_comm, mp_mat, mp_ierr); CHKERRA(mp_ierr)
        CALL MatSetSizes(mp_mat, m_nOwnedDofs, m_nOwnedDofs, PETSC_DECIDE, PETSC_DECIDE, mp_ierr); CHKERRA(mp_ierr)

        IF (m_numprocs > 1) THEN
            CALL MatSetType(mp_mat, MATMPIAIJ, mp_ierr); CHKERRA(mp_ierr)
            CALL MatMPIAIJSetPreallocation(mp_mat, m_maxFacesPerCell + 1, PETSC_NULL_INTEGER, m_maxFacesPerCell, PETSC_NULL_INTEGER, mp_ierr); CHKERRA(mp_ierr)
        ELSE
            CALL MatSetType(mp_mat, MATSEQAIJ, mp_ierr); CHKERRA(mp_ierr)
            CALL MatSeqAIJSetPreallocation(mp_mat, m_maxFacesPerCell + 1, PETSC_NULL_INTEGER, mp_ierr); CHKERRA(mp_ierr)
        ENDIF
    END SUBROUTINE allocate_petsc_mat
    

    !...set PETSc global matrix and rhs vector (matrix-assembled method)
    SUBROUTINE set_petsc_matrix_n_rhs(A, rhs)
        IMPLICIT NONE
        Mat A
        Vec rhs

        !...local variables
        INTEGER(int_p) :: i, j, mat, f, currf, ibface, bcon, pol, sout, mat1
        INTEGER(int_p) :: ic_1, ic_2, pol1, pol2

        REAL(real_p) :: t_12, r_12, t_21, r_21, jface1, jface2, jface11, jface22
        REAL(real_p) :: velener_de, velener_r, n_omega1, n_omega2, disfun1

        REAL(real_p) :: temper, vol, gamma, sdotn, insdotn, vomega, vel, alpha, &
            rmuval, rxival, retval, inrmuval, inrxival, inretval, tempwall, &
            actin, actout, beta, mfp, diag_val, offDiag_val, jout, temper1

        PetscScalar, DIMENSION(m_maxFacesPerCell + 1) :: mat_values
        PetscInt, DIMENSION(m_maxFacesPerCell + 1) :: mat_colInds
        PetscInt rowIdx
        PetscScalar :: rhs_value, mat_value
        INTEGER(int_p) :: offDiag_count

        pol = polar_type(iband_ctx, p_ctx)
        vel = gpvel(iband_ctx, pol)
        vomega = omega(si_ctx)
        rmuval = rmu(si_ctx)
        rxival = rxi(si_ctx)
        retval = ret(si_ctx)
        inrmuval = inrmu(si_ctx)
        inrxival = inrxi(si_ctx)
        inretval = inret(si_ctx)

        !...initialize matrix and rhs
        CALL VecZeroEntries(rhs, mp_ierr)
        CALL MatZeroEntries(A, mp_ierr)

        !... main loop finds the coefficients of the matrix A, and the rhs b
        DO i = 1, m_nOwnedDofs
            mat = cell_type(i)          !! material_type
            !IF (mat == ALUMINUM) CYCLE
            temper = tnot(i)

            !f = checkpolar(mat,iband_ctx, p_ctx)
            !IF (f == 0) GO TO 100

            vol = volcell(i) * vomega

            CALL relaxtime(pol, mat, iband_ctx, temper, beta)
            mfp = vel/beta
            
            rowIdx = lc_to_pgc(i) - 1 !... row index = global id of cell i (note PETSc uses 0-based index for both C and Fortran)
            mat_values(:) = 0.0
            rhs_value = 0.0
            offDiag_count = 0 !...number of off-diagonal terms of row i to be assembled

            DO j = 1,nfcell(i)
                currf = lcf(i,j) !... face id of face j of cell i
                gamma = mfp * areaf(currf)/vol

                IF (threed) THEN
                    sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival + vecfz(currf)*retval        
                    insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival+ vecfz(currf)*inretval
                ELSE 
                    sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival     
                    insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival
                ENDIF
                ! sdotn = sdotn_array(currf, si_ctx)
                ! insdotn = insdotn_array(currf, si_ctx)
                
                IF (bface(currf) == 0) THEN
                    alpha = normdir(i,j)
                    sdotn = alpha * sdotn
                    insdotn = alpha * insdotn
                    
                    !... column indices for diagonal term and off-diagonal term associated with face j
                    mat_colInds(1) = rowIdx
                    mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                    offDiag_count = offDiag_count + 1
                    mat_colInds(1 + offDiag_count) = lc_to_pgc(lcc(j,i)) - 1
                    
                    mat_values(1 + offDiag_count) = - Max(0.0,-sdotn)*gamma*((insdotn/(sdotn+tiny)))
                ELSE
                    !...Boundary face
                    ibface = f_to_bf(currf)
                    bcon = bctype(ibface)
                    tempwall = temp_bc(ibface)
                    
                    CALL iwalls(pol,mat, iband_ctx,tempwall,actin)
                    
                    jout = zero 
                    
                    SELECT CASE(bcon)
                    CASE(ADIA)
                        jout = jfacepos(ibface,pol,iband_ctx)

                        !... find the reflection direction
                        CALL reflec(si_ctx,currf,sdotn,sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)

                        !... assemble PETSc matrix
                        mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + dos*intensity(i,sout,icount_ctx,pol) * gamma * &
                                Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) &
                                +((one-dos)*jout) * gamma * Max(0.0,-sdotn)* &
                                ((insdotn/(sdotn+tiny)))/pi

                    CASE(ISOT)
                        !... assemble PETSc matrix
                        mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    CASE (HFLUX) !!! assuming a virtual block of aluminium!!
                        !... assemble PETSc matrix
                        mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    CASE (ALSI)
                        !... assemble PETSc matrix
                        mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))
                    
                    !!! symmetry boundary
                    CASE(SYMM)
                        jout = jfacepos(ibface,pol,iband_ctx)
                        CALL reflec(si_ctx,currf,sdotn,sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)

                        !... assemble PETSc matrix
                        mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + intensity(i,sout,icount_ctx,pol) * gamma * &
                                        Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    !!! INTERFACE....     
                    ! CASE (CONJ)
                    !     ic_1 = lfc(currf, 1)
                    !     ic_2 = lfc(currf, 2)
                        
                    !     pol1 = pol
                    !     pol2 =pol
                        
                    !     IF (cell_type(ic_1) /= cell_type(ic_2)) THEN
                    !         pol1 = mat_polar(cell_type(ic_1),pol)
                    !         pol2 = mat_polar(cell_type(ic_2),pol)
                    !     ENDIF
                    
                    !     CALL trans_coeff (iband_ctx,pol1,pol2, ic_1, ic_2,ibface, t_12)  

                    !     !!! file name changed, compilation done according to need              
                        
                    !     !t_12 = 0.5
                    !     r_12 = 1-t_12
                    !     t_21 = 1 - t_12
                    !     !t_21 = t_12
                    !     r_21 = 1 -t_21
                                                
                    !     jface1 = jfacepos(ibface,pol,iband_ctx)
                    !     jface2 = jfaceneg(ibface,pol,iband_ctx)
                    !     jface11 = jfacepos(ibface,pol1,iband_ctx)
                    !     jface22 = jfaceneg(ibface,pol2,iband_ctx)

                    !     IF (centfreq(iband_ctx) .lt. wmaxtas) THEN        
                    !         IF (cell_type(ic_1) == SILICON) THEN
                    !             jface11 = jfacepos(ibface,LA_1,iband_ctx) + jfacepos(ibface, TA_1,iband_ctx)
                    !             jface22 = jfaceneg(ibface,LA_2,iband_ctx) + jfaceneg(ibface, TA_2,iband_ctx)
                    !         ELSE
                    !             jface11 = jfacepos(ibface,LA_2,iband_ctx) + jfacepos(ibface, TA_2,iband_ctx)
                    !             jface22 = jfaceneg(ibface, LA_1,iband_ctx) + jfaceneg(ibface,TA_1,iband_ctx)
                            
                    !         ENDIF
                    !         jface1 = jface11
                    !         jface2 = jface22
                        
                    !         disfun1= one/(exp(hobol*centfreq(iband_ctx)/(tempwall)) - one)
                    
                    !         IF (cell_type(i) == SILICON) THEN
                    !             n_omega1 = (wavenumber(iband_ctx, LA_1))**2/(2* pi**2*(gpvel(iband_ctx,LA_1)))
                    !             n_omega2 =two*(wavenumber(iband_ctx, TA_1))**2/(2* pi**2*(gpvel(iband_ctx,TA_1))) 
                    !             velener_de = gpvel(iband_ctx,LA_1)*n_omega1*disfun1 + gpvel(iband_ctx,TA_1)*n_omega2*disfun1  
                    !         ELSE
                    !             n_omega1 = (wavenumber(iband_ctx, LA_2))**2/(2* pi**2*(gpvel(iband_ctx,LA_2)))
                    !             n_omega2 = two*(wavenumber(iband_ctx, TA_2))**2/(tiny+2* pi**2*(gpvel(iband_ctx,TA_2))) 
                    !             velener_de = gpvel(iband_ctx,LA_2)*n_omega1*disfun1 + gpvel(iband_ctx,TA_2)*n_omega2*disfun1
                    !         ENDIF
                    !         n_omega1 = (wavenumber(iband_ctx, pol))**2/(2* pi**2*(gpvel(iband_ctx,pol)))
                    !         IF ((pol .eq. TA_2) .or. (pol .eq. TA_1))THEN
                    !             n_omega1 = two*n_omega1 
                    !         ENDIF 
                    !         velener_r = (gpvel(iband_ctx,pol)*n_omega1*disfun1)/velener_de
                    !     ELSE
                    !         velener_r =  1.0
                    !     ENDIF  
                    !     !!!!!!! swartz modification ends

                    !     IF (ic_1 == i) THEN
                    !         !... assemble PETSc matrix
                    !         mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                    !         !... assemble PETSc rhs vector
                    !         rhs_value = rhs_value + velener_r*((jface1*r_12)/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) + &
                    !                     (jface22*t_21/pi)*gamma*Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))*velener_r

                    !     ELSE
                    !         !... assemble PETSc matrix
                    !         mat_values(1) = mat_values(1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                    !         !... assemble PETSc rhs vector
                    !         rhs_value = rhs_value + velener_r*((jface2*r_21)/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) + &
                    !                     (jface11*t_12/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))*velener_r
                    !     ENDIF                                                                  
                    END SELECT !! case of bcon (type of BCs of face j of cell i)
                ENDIF !Bounadry or Interior Face  
                
                !  scs_tran(i,pol) = scs_tran(i,pol)*gamma       !!! added 28th aug 
            ENDDO ! Faces of cell    loop j
          
            !...set to PETSc matrix and rhs vector
            CALL MatSetValues(A, 1, rowIdx, 1+offDiag_count, mat_colInds, mat_values, ADD_VALUES, mp_ierr)
            CALL VecSetValue(rhs, rowIdx, rhs_value, ADD_VALUES, mp_ierr)

            !!!!!13thaug$$$!!!!!!!!!
        100 ENDDO       !! 1st cell loop

        !... now compute the contribution of the first and the third terms of eq (3.17)
        DO i = 1,m_nOwnedDofs
            mat1 = cell_type(i)
            ! IF  (mat1== ALUMINUM) CYCLE 
            temper1= tnot(i)
            !scs(i) = scs(i) + scs_tran(i,pol)       !!! added 28th aug 
            !f = checkpolar(mat1,iband_ctx, p_ctx)
            !IF (f /= 0) THEN
                CALL relaxtime(pol, mat1, iband_ctx, temper1, beta)
                !beta = 1.0D-6
                CALL iwalls(pol,mat1, iband_ctx,temper1,actout)

                !... assemble PETSc matrix
                rowIdx = lc_to_pgc(i) - 1 !rowIdx = lc_to_gc(i) - 1
                mat_value = (idt/beta + one)
                CALL MatSetValue(A, rowIdx, rowIdx, mat_value, ADD_VALUES, mp_ierr);

                !... assemble PETSc rhs vector
                rhs_value = actout + idt *intensityone(i,si_ctx,icount_ctx,pol)/beta
                CALL VecSetValue(rhs, rowIdx, rhs_value, ADD_VALUES, mp_ierr)
                    
            !ENDIF
        ENDDO

        !...assemble matrix and rhs
        CALL MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, mp_ierr)
        CALL MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, mp_ierr)
        CALL VecAssemblyBegin(rhs, mp_ierr)
        CALL VecAssemblyEnd(rhs, mp_ierr)
    END SUBROUTINE set_petsc_matrix_n_rhs


    !...set s = -intensity of previous iteration for option 2 of petsc
    SUBROUTINE set_sol_prev_itsy_petsc(s)
        IMPLICIT NONE
        Vec s

        !...internal variables
        INTEGER(int_p) :: i
        PetscScalar, Pointer, DIMENSION(:) :: ss
        PetscErrorCode p_err

        !...point to petsc vector allowing ss to read/write to s
        CALL VecGetArrayF90(s, ss, p_err); CHKERRA(p_err)
        DO i = 1,m_nOwnedDofs
            ss(i) = -intensity(i, si_ctx, icount_ctx, p_ctx)
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(s, ss, p_err); CHKERRA(p_err)
        
    END SUBROUTINE set_sol_prev_itsy_petsc
END MODULE petsc_matAssembled
