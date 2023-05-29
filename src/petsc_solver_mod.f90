! \brief parent module coverings variables and methods related to PETSc
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! Setup PETSc variables and methods for solving the linear systems of intensity.
! This module is used by both methods of matrix assembled, assembly free (cell-based, combined
! band+cell-based, and batched cell-based parallelizations)
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
MODULE petsc_solver
    USE PRECISIONS
    USE constants
    USE VARIABLES, ONLY: m_rank => rank_cell, &     !...rank in cell-based parallel
        m_numprocs => numprocs_cell, &              !...number of processes in cell-based parallel
        m_comm => comm_cell, &                      !...communicator in cell-based parallel
        tnot, omega, gpvel, rmu, rxi, ret, &
        inrmu, inrxi, inret, jfacepos, jfaceneg, centfreq, &
        polar_type, SILICON, GERMENIUM, ALUMINUM, LA_1, TA_1, &
        wavenumber, idt, dos, intensity, intensityone, &
        polar_type, np_max, ndir!, sdotn_array, insdotn_array
    USE GRID, ONLY: m_nOwnedDofs => ncells, &       !...number of owned dofs
        m_nLocalDofs => nCells_total, &             !...number of local dofs = nOwnedDofs + nGhostDofs
        m_nGhostDofs => nGhostCells, &              !...number of ghost dofs
        lc_to_gc, lc_to_pgc, &
        volcell, threed, vecfx, vecfy, vecfz, nfcell, lcf, areaf, bface, temp_bc, &
        bctype, f_to_bf, ADIA, ISOT, CONJ, HFLUX, ALSI, SYMM, lcc, normdir, &
        cell_type, lfc
    USE MPI
    USE PETSC
    USE omp_lib
    USE profiler

    IMPLICIT NONE

    !...variables for PETSc solvers
    Vec mp_sol              !...pointer to solution vector
    Vec mp_rhs              !...pointer to rhs vector
    Mat mp_mat              !...pointer to matrix (for both matrix-assembed or matrix-free approaches)
    KSP mp_ksp              !...ksp for solving system of linear equation associated with mp_mat and mp_rhs
    PC mp_pc                !...preconditioner
    PetscErrorCode mp_ierr  !...error code for PETSc subroutines

    !...constants for PETSc solvers
    KSPType, PARAMETER :: mp_KSPTYPE = KSPGMRES     !...ksp type
    PCType, PARAMETER :: mp_PCTYPE = PCBJACOBI      !...default preconditioner type (e.g. PCSOR, PCLU, PCJACOBI, PCNONE)
    REAL(real_p), PARAMETER :: mp_RTOL = 1.0D-9     !...relative tolerance
    REAL(real_p), PARAMETER :: mp_ATOL = 1.0D-15    !...absolute tolerance
    INTEGER(int_p), PARAMETER:: mp_MAXITER = 3000   !...max number of iterations
    INTEGER(int_p) :: m_maxFacesPerCell             !...max number of faces per cell, used to allocal m_matOffDiag and array to set mp_matPC

    PetscInt :: m_nIts_total, m_nIts

    REAL(real_p) :: aMatMatvec_max, copyMatMult_max

    LOGICAL :: m_matPC_flag !...flag if PETSc matrix used for block-jacobi / SOR preconditioner was allocated
    Mat mp_matPC            !...pointer to PETSc matrix for block-jacobi / SOR precondition

    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_sorX       !...current vector used in SOR preconditioner
    
    INTEGER(int_p) :: iband_ctx, si_ctx, p_ctx, icount_ctx  !...get from disord
    INTEGER(int_p) :: iter_ctx
    
    CONTAINS
    !...determine max number faces per cell
    SUBROUTINE determine_max_faces_per_cell(max_nfaces)
        IMPLICIT NONE

        INTEGER(int_p), INTENT(OUT) :: max_nfaces
        INTEGER(int_p) :: i

        !...loop over cells (local) to find max faces per cell (this can be used for d_nz in matrix allocation)
        max_nfaces = 0
        DO i = 1,m_nOwnedDofs
            IF (nfcell(i) > max_nfaces) THEN
                max_nfaces = nfcell(i)
            ENDIF
        ENDDO
    END SUBROUTINE determine_max_faces_per_cell


    !...allocate, set type and size of the rhs vector, sol vector is a duplicate of rhs
    SUBROUTINE allocate_petsc_rhs_n_sol
        IMPLICIT NONE

        CALL VecCreate(m_comm, mp_rhs, mp_ierr); CHKERRA(mp_ierr)
        IF (m_numprocs > 1) THEN
            CALL VecSetType(mp_rhs, VECMPI, mp_ierr); CHKERRA(mp_ierr)
        ELSE
            CALL VecSetType(mp_rhs, VECSEQ, mp_ierr); CHKERRA(mp_ierr)
        ENDIF
        CALL VecSetSizes(mp_rhs, m_nOwnedDofs, PETSC_DECIDE, mp_ierr); CHKERRA(mp_ierr)
        CALL VecDuplicate(mp_rhs, mp_sol, mp_ierr); CHKERRA(mp_ierr)
    END SUBROUTINE allocate_petsc_rhs_n_sol


    !...create ksp, associate ksp with the matrix (for all 3 methods)
    SUBROUTINE create_petsc_ksp
        IMPLICIT NONE
        PetscReal :: rtol, abstol, dtol
        PetscInt :: maxits
        
        CALL KSPCreate(m_comm, mp_ksp, mp_ierr); CHKERRA(mp_ierr)
        CALL KSPSetType(mp_ksp, mp_KSPTYPE, mp_ierr); CHKERRA(mp_ierr) !...set default ksp is mp_KSPTYPE

        CALL KSPSetTolerances(mp_ksp, mp_RTOL, mp_ATOL, PETSC_DEFAULT_REAL, mp_MAXITER, mp_ierr); CHKERRA(mp_ierr)
        CALL KSPSetFromOptions(mp_ksp, mp_ierr); CHKERRA(mp_ierr) !...allow users to set ksp type at run time
        
        CALL KSPGetPC(mp_ksp, mp_pc, mp_ierr)
        CALL PCSetType(mp_pc, mp_PCTYPE, mp_ierr)
        CALL PCSetFromOptions(mp_pc, mp_ierr)

        !...echo tolerances
        CALL KSPGetTolerances(mp_ksp, rtol, abstol, dtol, maxits, mp_ierr)
        IF (m_rank .eq. 0) THEN
            WRITE(*,'(a,e15.3)') 'Inner relative tolerance = ', rtol
            WRITE(*,'(a,e15.3)') 'Inner absolute tolerance = ', abstol
            WRITE(*,'(a,i0)') 'Inner max number of iterations = ', maxits
        ENDIF

    END SUBROUTINE create_petsc_ksp


    !...deallocate PETSc vectors mp_rhs and mp_sol
    SUBROUTINE destroy_petsc_rhs_n_sol
        IMPLICIT NONE
        
        CALL VecDestroy(mp_rhs, mp_ierr); CHKERRA(mp_ierr)
        CALL VecDestroy(mp_sol, mp_ierr); CHKERRA(mp_ierr)
    END SUBROUTINE destroy_petsc_rhs_n_sol


    !...deallocate PETSc matrix
    SUBROUTINE destroy_petsc_mat
        IMPLICIT NONE
        
        CALL MatDestroy(mp_mat, mp_ierr); CHKERRA(mp_ierr)
    END SUBROUTINE destroy_petsc_mat


    !...deallocate PETSc ksp
    SUBROUTINE destroy_petsc_ksp
        IMPLICIT NONE

        CALL KSPDestroy(mp_ksp, mp_ierr); CHKERRA(mp_ierr)
    END SUBROUTINE destroy_petsc_ksp

END MODULE petsc_solver
