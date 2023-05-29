!...this subroutine solve for intensity using batched amat assembly-free approach
#include <petsc/finclude/petscksp.h>
SUBROUTINE disord_aMat_bch
    USE precisions
    USE constants
    USE grid, Only : ncells, vecfx, vecfy, vecfz, lfc, temp_bc, bctype, &
        CONJ, nbcfaces, bf_to_f, cell_type, threed, lc_to_gc, lc_to_pgc
    USE variables
    USE amat_batched

    IMPLICIT NONE

    !... pointer to Petsc solution vector
    PetscScalar, Pointer, DIMENSION(:) :: sol_v

    INTEGER(int_p) :: i, j, currf, ic, ifcb, p
    INTEGER(int_p) :: icell
    INTEGER(int_p) :: si, iband, pol, mat, f
    INTEGER(int_p) :: idxOffset

    INTEGER(int_p) :: icount  !!!  intensity band index

    REAL(real_p) :: sdotn, insdotn
    REAL(real_p) :: vel
    REAL(real_p) :: actin, tempnum, beta, jout, ftemp, jin , temp_flux

    !! interface variables
    INTEGER (int_p) :: ic_1, ic_2, pol1, pol2
    REAL(real_p) :: t_12, r_12, t_21, r_21, jface1, jface2

    INTEGER(int_p) :: b
    
    INTEGER(int_p) :: count
    REAL(real_p) :: normRhs

    IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting bsnord_disord"

    DO iband = lowband, highband
        icount = (iband - lowband + 1)
        
        !...context data of band (same for all directions and polarizations)
        iband_ctx = iband
        icount_ctx = icount

        !...loop over number of batches to solve the "batched system" (i.e. solve all non-zero systems in batch i at once)
        DO i = 1, nBatchesPerBand(iband)
            !...get the context data of each non-zero system in batch i
            DO j = 1,m_batchCount
                si_ctx_bch(j) = si_nonzeroSystemPerBand(iband, (i-1)*m_batchCount + j)
                p_ctx_bch(j) = p_nonzeroSystemPerBand(iband, (i-1)*m_batchCount + j)
            ENDDO
            
            !...compute "elemental" matrices and rhs vector
            !CALL profiler_start(matRhsComp_timer)
            CALL compute_aMat_matrix_n_rhs_bch
            CALL build_petsc_diagonal_block_matrix_bch
            !CALL profiler_stop(matRhsComp_timer)

            !CALL profiler_start(preSolve_timer)
            CALL modify_aMat_rhs_for_dx_solving_bch_1   !...2022.03.19: also compute the norm of rhs of individual system in the batch

            !...set m_rhs (owned by aMat) to PETSc vector mp_rhs (this includes VecAssemblyBegin/End)
            CALL set_petsc_rhs_bch(mp_rhs)
            !CALL profiler_stop(preSolve_timer)
            
            !...solve using amat approach
            CALL KSPSetOperators(mp_ksp, mp_mat, mp_mat, mp_ierr) !...re-compute preconditioner (if any)
            CALL KSPSolve(mp_ksp, mp_rhs, mp_sol, mp_ierr)
            
            !... get pointer to solution vector; NOTE: sol_v use 1-based index as regular Fortran array
            CALL VecGetArrayF90(mp_sol, sol_v, mp_ierr); CHKERRA(mp_ierr)

            !... put back the solution of "batched system" in intensity array
            DO j = 1,m_batchCount
                si = si_ctx_bch(j)
                p = p_ctx_bch(j)
                DO icell = 1,ncells
                    !IF(cell_type(icell) == ALUMINUM) CYCLE
                    intensity(icell, si, icount, p) = intensity(icell, si, icount, p) + sol_v(m_ocSid2OcMap(icell,j))
                ENDDO
            ENDDO
        ENDDO !...i = 1,nBatchesPerBand(iband)
    ENDDO ! band loop ends. more stream lining, so that flux and gc is calculated once in iteration

    !... so far, we solved for all intensity, then just integrating the intensity to find the flux
    ! Calculation of flux at boundary faces
    DO iband = lowband, highband 
        icount = (iband - lowband + 1)
        DO p = 1, np_max      !! loop over polarization
            pol = polar_type(iband, p)
            IF (pol .eq. 0) CYCLE
            vel = gpvel(iband, pol) !! velocity depends on band, polarization

            !...loop over boundary faces
            DO ifcb = 1, nbcfaces
                !...Han? why not using bf_to_c
                currf = bf_to_f(ifcb)
                ic = lfc(currf,1)
                ! IF (bctype(ifcb) == ADIA .and. celL_type(ic) == ALUMINUM) CYCLE
                IF (bctype(ifcb) /= CONJ) THEN
                !IF(bctype(ifcb) == ALSI .and. cell_type(ic)== ALUMINUM) ic = lfc(currf,2)  
                    !IF (cell_type(ic)== ALUMINUM) ic = lfc(currf,2)     
                    ftemp = temp_bc(ifcb)
                    
                    mat = cell_type(ic)
                    temp_flux= zero
                    !f = checkpolar(mat,iband, p)
                    !IF (f /= 0) THEN 
                        CALL iwalls(pol,mat, iband,ftemp,actin)
                        jout = pi*actin
                        jin = zero
                        DO si = 1,ndir            
                            !sdotn =  vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)  + vecfz(currf) * ret(si)
                            !insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si) +  vecfz(currf) * inret(si)
                            IF (threed) THEN       
                                sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si) + vecfz(currf)*ret(si)        
                                insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)+ vecfz(currf)*inret(si)
                            ELSE 
                                sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)     
                                insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si)
                            ENDIF
                            jin = jin + max(0.0,sdotn)*intensity(ic,si,icount,pol)*(insdotn/sdotn)                               
                        ENDDO   !! dir
                        flux(ifcb,pol,iband) = (jin-jout)
                        flux_s(ifcb,pol,iband) = (jin-jout)
                        !flux(ifcb,pol,iband) = jfacepos(ifcb,pol,iband) - jout !(jin-jout)
                        !flux_s(ifcb,pol,iband) = jfacepos(ifcb,pol,iband) - jout !(jin-jout)
                    !ENDIF
                ELSEIF (bctype(ifcb) == CONJ) THEN
                    ! currf = bf_to_f(ifcb)
                    ! ic_1 = lfc(currf,1)
                    ! ic_2 = lfc (currf,2)
                    ! pol1 = pol
                    ! pol2 =pol
                    ! IF (cell_type(ic_1) /= cell_type(ic_2)) THEN
                    !     pol1 = mat_polar(cell_type(ic_1),pol)
                    !     pol2 = mat_polar(cell_type(ic_2),pol)
                    ! ENDIF
                                
                    ! ! CALL trans_coeff_dmm1 (iband,pol1,pol2,ic_1, ic_2,ifcb, t_12) 
                    ! ! CALL trans_coeff_dmm (iband, pol1, ic_1, ic_2,ifbc, t_12)
                    ! !!! file name changed and compile according to need
                    ! CALL  trans_coeff (iband,pol1,pol2, ic_1, ic_2,ifcb, t_12)     
                        
                    ! r_12 = 1-t_12
                    ! t_21 = 1 - t_12
                    ! r_21 = 1 -t_21
                    ! jface2 = zero
                    ! jface1 = zero
                    ! DO si = 1,ndir            
                    !     sdotn = vecfx(currf) * rmu(si) + vecfy(currf) * rxi(si)  + &
                    !             vecfz(currf) * ret(si)
                    !     insdotn = vecfx(currf) * inrmu(si) + vecfy(currf) * inrxi(si) + &
                    !             vecfz(currf) * inret(si)
                    !     jface1 = jface1 + max(0.0,sdotn)*intensity(ic_1,si,icount,pol)*(insdotn/sdotn)
                    !     jface2 = jface2 - max(0.0,-sdotn)*intensity(ic_2,si,icount,pol)*(insdotn/sdotn)   
                    ! ENDDO   !! dir
                    ! flux(ifcb,pol,iband) = t_12*ABS(jface1) - t_21*ABS(jface2)   !! flux is the destination of gathered flux_s
                    ! flux_s(ifcb,pol,iband) = t_12*ABS(jface1) - t_21*ABS(jface2) !! flux_s is used to avoid some error of memory-access errro
                    CYCLE
                ENDIF 
            ENDDO       !! nbc faces
        ENDDO !...p = 1,np_max
                
        ! !! this is to find gc (integrated intensity over directions)
        ! !         write(*,*) "Starting integral intensity"
        ! !!!!!!         !!! integral intesity            
        DO ic = 1, ncells
            mat = cell_type(ic)
            !IF  (mat== ALUMINUM) CYCLE 
            DO p = 1, np_max
                pol = polar_type(iband,p)
                IF (pol .eq. 0) CYCLE
                !CALL relaxtime(pol, mat, iband, 100.0, beta)
                !f=  checkpolar(mat, iband, p)
                !IF (f /=0) THEN
                !                   
                tempnum = zero
                DO si = 1, ndir
                    tempnum = tempnum + intensity(ic,si,icount,pol)*omega(si)
                ENDDO ! again direction

                gc(ic,pol,iband) = tempnum
                
                !gc(ic,iband) = gc(ic,iband) + tempnum
                !ENDIF
            ENDDO       !! polarization loop
            !gc(ic,iband) = tempnum
        ENDDO  ! loop ic, over ncells
        !!!!!!!!!!!
        !WRITE (*,*) "band done", iband
        !END DO                   !! for polarization
    ENDDO ! Bands iband

    IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing bsnord disord"
    10 format(4(1x,E14.6))
END SUBROUTINE disord_aMat_bch