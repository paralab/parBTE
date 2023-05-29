! \brief module to solve using assembly-free aMat, with batching (i.e. solve multiple linear systems simultaneously)
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! Solve the linear systems of intensity using assembly-free method, with batching approach
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
MODULE amat_batched
    USE petsc_solver
    USE amat_parent

    IMPLICIT NONE

    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_ui, m_vi           !...to get intermediate solution in batch method

    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: si_ctx_bch, p_ctx_bch
    INTEGER(int_p) :: m_batchCount                                  !...the number of systems Ax = b to be solved
    INTEGER(int_p) :: m_nOwnedDofs_bch                              !...number of owned dofs
    INTEGER(int_p) :: m_nLocalDofs_bch                              !...number of local (i.e. ghost included) dofs
    INTEGER(int_p) :: m_totalSend_bch                               !...total values to send to other processes
    INTEGER(int_p) :: m_totalRecv_bch                               !...total values to be received (= number of ghost)
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendCellIds_bch  !..."spread" of cell ids for batched sending
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendCounts_bch   !...counts of sending in batched communication
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendOffset_bch   !...offset of sendCounts
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_recvCounts_bch   !...counts of receiving in batched communication
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_recvOffset_bch   !...offset of recvCounts
    INTEGER(int_p), DIMENSION(:,:), ALLOCATABLE :: m_ocSid2GcMap    !...map from owned cell id and system id to global cell id of batched system
    INTEGER(int_p), DIMENSION(:,:), ALLOCATABLE :: m_lcSid2LcMap  
    INTEGER(int_p), DIMENSION(:,:), ALLOCATABLE :: m_ocSid2OcMap  
    
    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_rhsNorm_bch, m_rhsNormSum_bch  !...norm of rhs for each system in the batch (local, reduced sum)
    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_resNorm_bch, m_resNormSum_bch  !...norm of residual for each system in the batch (local, reduced sum)

    CONTAINS
    
    !...batched version of init_aMat
    SUBROUTINE init_aMat_bch(nSystemsPerBatch)
        IMPLICIT NONE
        INTEGER(int_p), INTENT(in) :: nSystemsPerBatch

        INTEGER(int_p) :: i, ierr

        !...get the batchCount from external variables
        m_batchCount = nSystemsPerBatch

        !...set the batch sizes of corresponding number owned dofs and local dofs
        m_nOwnedDofs_bch = m_batchCount * m_nOwnedDofs
        m_nLocalDofs_bch = m_batchCount * m_nLocalDofs

        !...allocate space for vectors u and v including ghost cells
        ALLOCATE(m_u(m_nLocalDofs_bch), stat=ierr)
        ALLOCATE(m_v(m_nLocalDofs_bch), stat=ierr)

        ALLOCATE(m_ui(m_nLocalDofs_bch), stat=ierr)
        ALLOCATE(m_vi(m_nLocalDofs_bch), stat=ierr)

        !...diagonal terms of rows owned by me
        ALLOCATE(m_matDiag(m_nOwnedDofs_bch), stat=ierr)
        
        !...determine max number of faces per cell
        CALL determine_max_faces_per_cell(m_maxFacesPerCell)

        !...off-diagonal terms of rows owned by me
        ALLOCATE(m_matOffDiag(m_maxFacesPerCell, m_nOwnedDofs_bch), stat=ierr) !...Fortran uses column-major order (consecutive elements of a column are stored next to each other)

        !...rhs vector and array of corresponding global indices
        ALLOCATE(m_rhs(m_nOwnedDofs_bch), stat=ierr)
        ALLOCATE(m_idxs(m_nOwnedDofs_bch), stat=ierr)

        !...allocate space for variables holding context
        AlLOCATE(si_ctx_bch(m_batchCount), stat=ierr)
        AlLOCATE(p_ctx_bch(m_batchCount), stat=ierr)

        !...allocate map from owned cell id and system id to local and global cell id
        ALLOCATE(m_ocSid2GcMap(m_nOwnedDofs, m_batchCount), stat=ierr)
        ALLOCATE(m_lcSid2LcMap(m_nLocalDofs, m_batchCount), stat=ierr)
        ALLOCATE(m_ocSid2OcMap(m_nOwnedDofs, m_batchCount), stat=ierr)

        !...allocate array storing norm of rhs of each system in the batch
        ALLOCATE(m_rhsNorm_bch(m_batchCount), stat=ierr)
        ALLOCATE(m_rhsNormSum_bch(m_batchCount), stat=ierr)

        !...allocate array storing norm of residual of each system in the batch
        ALLOCATE(m_resNorm_bch(m_batchCount), stat=ierr)
        ALLOCATE(m_resNormSum_bch(m_batchCount), stat=ierr)

        !...determine dep/ind cells
        CALL identifyIndCells

        !...profiling
        !CALL profiler_init(m_mvTimer)
        !CALL profiler_init(m_precondTimer)

        !...flag if block-jacobi matrix is allocated
        m_matPC_flag = .FALSE.
        
    END SUBROUTINE init_aMat_bch

    !...batched version of finalize_aMat
    SUBROUTINE finalize_aMat_bch
        IMPLICIT NONE

        IF (ALLOCATED(m_sendCellIds)) DEALLOCATE(m_sendCellIds)
        IF (ALLOCATED(m_sendCellCounts)) DEALLOCATE(m_sendCellCounts)
        IF (ALLOCATED(m_sendCellOffset)) DEALLOCATE(m_sendCellOffset)
        IF (ALLOCATED(m_recvCellCounts)) DEALLOCATE(m_recvCellCounts)
        IF (ALLOCATED(m_recvCellOffset)) DEALLOCATE(m_recvCellOffset)

        IF (ALLOCATED(m_sendBuf)) DEALLOCATE(m_sendBuf)
        IF (ALLOCATED(m_recvBuf)) DEALLOCATE(m_recvBuf)
        IF (ALLOCATED(m_sendRankIds)) DEALLOCATE(m_sendRankIds)
        IF (ALLOCATED(m_recvRankIds)) DEALLOCATE(m_recvRankIds)
        IF (ALLOCATED(m_request)) DEALLOCATE(m_request)
        IF (ALLOCATED(m_u)) DEALLOCATE(m_u)
        IF (ALLOCATED(m_v)) DEALLOCATE(m_v)
        IF (ALLOCATED(m_ui)) DEALLOCATE(m_ui)
        IF (ALLOCATED(m_vi)) DEALLOCATE(m_vi)
        IF (ALLOCATED(m_indCells)) DEALLOCATE(m_indCells)
        IF (ALLOCATED(m_depCells)) DEALLOCATE(m_depCells)
        IF (ALLOCATED(m_matDiag)) DEALLOCATE(m_matDiag)
        IF (ALLOCATED(m_matOffDiag)) DEALLOCATE(m_matOffDiag)
        IF (ALLOCATED(m_rhs)) DEALLOCATE(m_rhs)
        IF (ALLOCATED(m_idxs)) DEALLOCATE(m_idxs)

        IF (ALLOCATED(m_sendCellIds_bch)) DEALLOCATE(m_sendCellIds_bch)
        IF (ALLOCATED(m_sendCounts_bch)) DEALLOCATE(m_sendCounts_bch)
        IF (ALLOCATED(m_sendOffset_bch)) DEALLOCATE(m_sendOffset_bch)
        IF (ALLOCATED(m_recvCounts_bch)) DEALLOCATE(m_recvCounts_bch)
        IF (ALLOCATED(m_recvOffset_bch)) DEALLOCATE(m_recvOffset_bch)

        IF (ALLOCATED(p_ctx_bch)) DEALLOCATE(p_ctx_bch)
        IF (ALLOCATED(si_ctx_bch)) DEALLOCATE(si_ctx_bch)

        IF (ALLOCATED(m_ocSid2GcMap)) DEALLOCATE(m_ocSid2GcMap)
        IF (ALLOCATED(m_lcSid2LcMap)) DEALLOCATE(m_lcSid2LcMap)
        IF (ALLOCATED(m_ocSid2OcMap)) DEALLOCATE(m_ocSid2OcMap)

        IF (ALLOCATED(m_rhsNorm_bch)) DEALLOCATE(m_rhsNorm_bch)
        IF (ALLOCATED(m_rhsNormSum_bch)) DEALLOCATE(m_rhsNormSum_bch)

        IF (ALLOCATED(m_resNorm_bch)) DEALLOCATE(m_resNorm_bch)
        IF (ALLOCATED(m_resNormSum_bch)) DEALLOCATE(m_resNormSum_bch)

        IF (m_matPC_flag) THEN
            CALL MatDestroy(mp_matPC, mp_ierr); CHKERRA(mp_ierr)
        ENDIF
        CALL MatDestroy(mp_mat, mp_ierr) !...although we create mp_mat as a shell, we should destroy it
    END SUBROUTINE finalize_aMat_bch

    !...batched version of allocate_amat_rhs_n_sol
    SUBROUTINE allocate_amat_rhs_n_sol_bch
        IMPLICIT NONE
        
        CALL VecCreate(m_comm, mp_rhs, mp_ierr); CHKERRA(mp_ierr)
        IF (m_numprocs > 1) THEN
            CALL VecSetType(mp_rhs, VECMPI, mp_ierr); CHKERRA(mp_ierr)
        ELSE
            CALL VecSetType(mp_rhs, VECSEQ, mp_ierr); CHKERRA(mp_ierr)
        ENDIF

        CALL VecSetSizes(mp_rhs, m_nOwnedDofs_bch, PETSC_DECIDE, mp_ierr); CHKERRA(mp_ierr)
        CALL VecSetFromOptions(mp_rhs, mp_ierr); CHKERRA(mp_ierr)

        CALL VecDuplicate(mp_rhs, mp_sol, mp_ierr); CHKERRA(mp_ierr)

    END SUBROUTINE allocate_amat_rhs_n_sol_bch

    !...batched version of allocate_petsc_matshell
    SUBROUTINE allocate_petsc_matshell_bch
        IMPLICIT NONE

        !...create matrix shell with batch size
        CALL MatCreateShell(m_comm, m_nOwnedDofs_bch, m_nOwnedDofs_bch, PETSC_DECIDE, PETSC_DECIDE, &
            PETSC_NULL_INTEGER, mp_mat, mp_ierr)

        !...define matrix-vector multiplication operator
        CALL MatShellSetOperation(mp_mat, MATOP_MULT, MatMult_amat_bch, mp_ierr)

        !...define get matrix diagonal operator, used in pc jacobi
        CALL MatShellSetOperation(mp_mat, MATOP_GET_DIAGONAL, MatGetDiagonal_amat_bch, mp_ierr)

        !...define get block diagonal operator, used in pc bjacobi
        CALL MatShellSetOperation(mp_mat, MATOP_GET_DIAGONAL_BLOCK, MatGetDiagonalBlock_wrapper_bch, mp_ierr)

        !...define get sor sweep, used in pc sor
        CALL MatShellSetOperation(mp_mat, MATOP_SOR, MatSOR_wrapper_bch, mp_ierr)

    END SUBROUTINE allocate_petsc_matshell_bch

    !...create ksp, associate ksp with the matrix (for all 3 methods)
    SUBROUTINE create_petsc_ksp_bch
        IMPLICIT NONE
        PetscReal :: rtol, abstol, dtol
        PetscInt :: maxits

        CALL KSPCreate(m_comm, mp_ksp, mp_ierr); CHKERRA(mp_ierr)
        CALL KSPSetType(mp_ksp, mp_KSPTYPE, mp_ierr); CHKERRA(mp_ierr) !...set default ksp is mp_KSPTYPE

        !...use user-defined convergence test
        CALL KSPSetConvergenceTest(mp_ksp, MyKSPConverged, 0, PETSC_NULL_FUNCTION, mp_ierr)

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

    END SUBROUTINE create_petsc_ksp_bch

    !...batched version of buildScatterMap
    SUBROUTINE buildScatterMap_bch
        IMPLICIT NONE

        INTEGER(int_p) :: i, j, counter1, counter2, ierr
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: localCellCounts, localCellOffset !...number of cells owned by each rank
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: ghostCounts, ghostOffset !...number of ghost cells owned by other rank
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: ghostGIds !...list of global Ids of ghost cells
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: ghostOwner !...owner of ghost cells
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: recvCounts, recvOffset !...number of owned cells that I need to send to each rank
        INTEGER(int_p) :: nRecvs
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: recvBuff

        INTEGER(int_p) :: globalId, localId
        INTEGER(int_p), PARAMETER :: root_proc = 0
        INTEGER :: mpierr
        INTEGER(int_p) :: r, b, idxOffset

        !...variables hold number of owned cells across ranks
        ALLOCATE(localCellCounts(0:m_numprocs-1), stat=ierr)
        ALLOCATE(localCellOffset(0:m_numprocs-1), stat=ierr)

        ALLOCATE(ghostCounts(0:m_numprocs-1), stat=ierr)
        ALLOCATE(ghostOffset(0:m_numprocs-1), stat=ierr)

        ALLOCATE(recvCounts(0:m_numprocs-1), stat=ierr)
        ALLOCATE(recvOffset(0:m_numprocs-1), stat=ierr)

        !...communicate to get m_nOwnedDofs of all other ranks
        CALL MPI_ALLGATHER(m_nOwnedDofs, 1, MPI_INTEGER, localCellCounts, 1, MPI_INTEGER, m_comm, mpierr)

        !...compute localCellOffset
        localCellOffset(0) = 0
        DO i = 1,m_numprocs-1
            localCellOffset(i) = localCellOffset(i-1) + localCellCounts(i-1)
        ENDDO

        !...start and end of global cell id (in pgc) owned by my rank
        m_globalDofBegin = localCellOffset(m_rank) + 1 !...pgc uses 1-based index
        m_globalDofEnd = localCellOffset(m_rank) + localCellCounts(m_rank) !...pgc uses 1-based index
        
        !...list of ghost ids in terms of "pgc"
        ALLOCATE(ghostGIds(m_nGhostDofs), stat=ierr)
        DO i = 1,m_nGhostDofs
            ghostGIds(i) = lc_to_pgc(m_nOwnedDofs + i)
        ENDDO

        !...determine the owner of ghost cells
        !...required: global cell ids are continuously partitioned across processors
        ALLOCATE(ghostOwner(m_nGhostDofs), stat=ierr)
        DO i = 1,m_nGhostDofs
            globalId = ghostGIds(i) !...global ID of ghost cell
            DO j = 0,m_numprocs-1
                IF ((globalId > localCellOffset(j)) .and. (globalId <= (localCellOffset(j) + localCellCounts(j)))) THEN
                    ghostOwner(i) = j !...rank id is 0-based index
                    EXIT
                ENDIF
            ENDDO
        ENDDO

        !...determine number of my ghost cells belong to each rank, e.g. I am rank 0, ghostCount(0) is always 0, 
        !... if ghostCounts(1) = 3 then I have 3 ghost cells whose owner is rank 1
        !... these are the number of intensity values that I will receive from each rank before MATVEC
        DO i = 0,m_numprocs-1
            ghostCounts(i) = 0
        ENDDO
        DO i = 1,m_nGhostDofs
            ghostCounts(ghostOwner(i)) = ghostCounts(ghostOwner(i)) + 1
        ENDDO
        ghostOffset(0) = 0
        DO i = 1,m_numprocs-1
            ghostOffset(i) = ghostOffset(i-1) + ghostCounts(i-1)
        ENDDO

        !...determine recvCounts (number of owned cells that I need to send to each rank before MATVEC)
        !...e.g. I am rank 0, if recvCounts(1) = 4 then I own 4 cells that rank 1 will need them as ghost cells
        CALL MPI_ALLTOALL(ghostCounts, 1, MPI_INTEGER, recvCounts, 1, MPI_INTEGER, m_comm, mpierr)
        
        recvOffset(0) = 0
        DO i = 1,m_numprocs-1
            recvOffset(i) = recvOffset(i-1) + recvCounts(i-1)
        ENDDO

        !...total number of cells that I need to send to other ranks who need them as ghost (before doing matvec)
        nRecvs = recvOffset(m_numprocs-1) + recvCounts(m_numprocs-1)

        !...total cells to be sent out
        m_totalSend_bch = nRecvs * m_batchCount !...total values sent in batch

        !...total cells to be received from
        m_totalRecv_bch = m_nGhostDofs * m_batchCount

        !...determine the cell Ids that I need to send to corresponding rank who need them as ghost
        !...first, allocate recv buffer to receive the global Ids (of owned cells that I need to send to other ranks before doing matvec)
        ALLOCATE(recvBuff(nRecvs), stat=ierr)

        !...send ghost Ids (global Ids) to ranks who own them so that the owner will know which cells that they need to send to corresponding rank before MATVEC
        CALL MPI_ALLTOALLV(ghostGIds, ghostCounts, ghostOffset, MPI_INTEGER, recvBuff, recvCounts, recvOffset, &
            MPI_INTEGER, m_comm, mpierr)

        !...put data in global variables representing scatter map (to be used in matvec subroutine)
        !...these variables are deallocated in finalize_aMat
        ALLOCATE(m_sendCellIds(nRecvs), stat=ierr)
        ALLOCATE(m_sendCellCounts(0:m_numprocs-1), stat=ierr)
        ALLOCATE(m_sendCellOffset(0:m_numprocs-1), stat=ierr)
        ALLOCATE(m_recvCellCounts(0:m_numprocs-1), stat=ierr)
        ALLOCATE(m_recvCellOffset(0:m_numprocs-1), stat=ierr)
        DO i = 0,m_numprocs-1
            m_sendCellCounts(i) = recvCounts(i) !...TODO: once the batched version works correctlly, we can get rid of this and compute directly sendCounts_batched here
            m_sendCellOffset(i) = recvOffset(i)
            m_recvCellCounts(i) = ghostCounts(i)
            m_recvCellOffset(i) = ghostOffset(i)
        ENDDO
        DO i = 1,nRecvs
            globalId = recvBuff(i) !...global ID of cell that I need to send to other ranks
            localId = globalId - localCellOffset(m_rank) !...convert to local cell id which is 1-based index
            m_sendCellIds(i) = localId
        ENDDO

        !...allocate buffers for send and recv in batched matvec
        IF (m_totalSend_bch > 0) THEN
            ALLOCATE(m_sendBuf(0:m_totalSend_bch-1), stat=ierr)
        ENDIF
        IF (m_totalRecv_bch > 0) THEN
            ALLOCATE(m_recvBuf(0:m_totalRecv_bch-1), stat=ierr)
        ENDIF

        !...identify ranks that I need to send to and ranks that I will receive from
        m_nProcsSend = 0
        m_nProcsRecv = 0
        DO i = 0,m_numprocs-1
            IF (m_sendCellCounts(i) > 0) m_nProcsSend = m_nProcsSend + 1
            IF (m_recvCellCounts(i) > 0) m_nProcsRecv = m_nProcsRecv + 1
        ENDDO
        ALLOCATE(m_sendRankIds(m_nProcsSend), stat=ierr)
        ALLOCATE(m_recvRankIds(m_nProcsRecv), stat=ierr)
        counter1 = 0
        counter2 = 0
        DO i = 0,m_numprocs-1
            IF (m_sendCellCounts(i) > 0) THEN
                counter1 = counter1 + 1
                m_sendRankIds(counter1) = i
            ENDIF
            IF (m_recvCellCounts(i) > 0) THEN
                counter2 = counter2 + 1
                m_recvRankIds(counter2) = i
            ENDIF
        ENDDO

        !...allocate space for storing requests (total n requests = nProcsSend + nProcsRecv)
        ALLOCATE(m_request(m_nProcsSend + m_nProcsRecv), stat=ierr)

        !...build maps from owned/local cell id and system id to global cell id, owned cell id and local cell id of batched system
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                m_ocSid2GcMap(i,b) = (m_batchCount * localCellOffset(m_rank)) + ((b - 1) * m_nOwnedDofs) + i
                m_ocSid2OcMap(i,b) = ((b - 1) * m_nOwnedDofs) + i
            ENDDO
            DO i = 1,m_nLocalDofs
                m_lcSid2LcMap(i, b) = (b - 1) * m_nLocalDofs + i
            ENDDO
        ENDDO

        !...following are for send/recv in batch
        ALLOCATE(m_sendCellIds_bch(nRecvs * m_batchCount))
        ALLOCATE(m_sendCounts_bch(0:m_numprocs-1), stat=ierr)
        ALLOCATE(m_sendOffset_bch(0:m_numprocs-1), stat=ierr)
        ALLOCATE(m_recvCounts_bch(0:m_numprocs-1), stat=ierr)
        ALLOCATE(m_recvOffset_bch(0:m_numprocs-1), stat=ierr)

        DO i = 0,m_numprocs-1
            m_sendCounts_bch(i) = m_sendCellCounts(i) * m_batchCount
            m_sendOffset_bch(i) = m_sendCellOffset(i) * m_batchCount
            m_recvCounts_bch(i) = m_recvCellCounts(i) * m_batchCount
            m_recvOffset_bch(i) = m_recvCellOffset(i) * m_batchCount
        ENDDO
        DO i = 1,m_nProcsSend 
            r = m_sendRankIds(i)
            DO b = 1,m_batchCount
                idxOffset = (b-1) * m_nLocalDofs
                DO j = 1,m_sendCellCounts(r)
                    m_sendCellIds_bch(m_sendOffset_bch(r) + (b-1)*m_sendCellCounts(r) + j) = &
                        m_sendCellIds(m_sendCellOffset(r) + j) + idxOffset
                ENDDO
            ENDDO
        ENDDO

        !...deallocate variables used in buildScatterMap
        IF (ALLOCATED(localCellCounts)) DEALLOCATE(localCellCounts)
        IF (ALLOCATED(localCellOffset)) DEALLOCATE(localCellOffset)
        IF (ALLOCATED(ghostCounts)) DEALLOCATE(ghostCounts)
        IF (ALLOCATED(ghostOffset)) DEALLOCATE(ghostOffset)
        IF (ALLOCATED(recvCounts)) DEALLOCATE(recvCounts)
        IF (ALLOCATED(recvOffset)) DEALLOCATE(recvOffset)
        IF (ALLOCATED(ghostGIds)) DEALLOCATE(ghostGIds)
        IF (ALLOCATED(ghostOwner)) DEALLOCATE(ghostOwner)
        IF (ALLOCATED(recvBuff)) DEALLOCATE(recvBuff)
    END SUBROUTINE buildScatterMap_bch

    !...batched version of ghost_receive_begin
    SUBROUTINE ghost_receive_begin_bch(vec)
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(IN) :: vec

        !...local variables
        INTEGER(int_p) :: i, j, r, counter, counter_1, req
        INTEGER :: mpierr
        INTEGER(int_p) :: batchId, idxOffset
        INTEGER(int_p), PARAMETER :: commTag = 1

        !...do nothing in sequential runs
        IF (m_numprocs == 1) RETURN

        !...send as a batch of values of v at (owned) cells to other ranks who need them as ghost
        counter = 0
        IF (m_totalSend_bch > 0) THEN
            !...put values of vec at positions indicated by m_sendCellIds_bch to sending buffer
            DO i = 1,m_totalSend_bch
                m_sendBuf(i-1) = vec(m_sendCellIds_bch(i))
            ENDDO
            !...send data to other processes
            DO i = 1,m_nProcsSend 
                r = m_sendRankIds(i) !...rank that I am going to send to
                CALL MPI_Isend(m_sendBuf(m_sendOffset_bch(r)), m_sendCounts_bch(r), MPI_REAL8, r, &
                    commTag, m_comm, req, mpierr)
                !... put req (generated by MPI_Isend) to the request list so that we can use later in MPI_Wait
                counter = counter + 1
                m_request(counter) = req
            ENDDO
        ENDIF

        !...receive values of v at ghost cells from their owners
        IF (m_totalRecv_bch > 0) THEN
            DO i = 1,m_nProcsRecv
                r = m_recvRankIds(i) !...rank that I received from
                CALL MPI_Irecv(m_recvBuf(m_recvOffset_bch(r)), m_recvCounts_bch(r), MPI_REAL8, r, &
                    commTag, m_comm, req, mpierr)
                !... put req (generated by MPI_Irecv) to the request list so that we can use later in MPI_Wait
                counter = counter + 1
                m_request(counter) = req
            ENDDO
        ENDIF
    END SUBROUTINE ghost_receive_begin_bch

    !...batched version of ghost_receive_end
    SUBROUTINE ghost_receive_end_bch(vec)
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(OUT) :: vec

        !...local variables
        INTEGER(int_p) :: i, j, k, b, count, idxOffset
        INTEGER :: mpierr

        IF (m_numprocs == 1) RETURN

        !...synchronize (wait) all requests of isend and irecv have started
        DO i = 1,(m_nProcsRecv + m_nProcsSend)
            CALL MPI_Wait(m_request(i), MPI_STATUS_IGNORE, mpierr)
        ENDDO
        
        !...put the received values to the ghost position of vec
        count = 0
        DO i = 1,m_nProcsRecv
            k = m_recvRankIds(i)
            DO b = 1,m_batchCount
                idxOffset = (b-1) * m_nLocalDofs
                DO j = 1,m_recvCellCounts(k)
                    count = count + 1
                    vec(idxOffset + m_recvCellOffset(k) + j + m_nOwnedDofs) = m_recvBuf(count-1)
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE ghost_receive_end_bch

    !...batched version of MatMult_amat
    SUBROUTINE MatMult_amat_bch(A, u, v, ierr)
        IMPLICIT NONE
        
        Mat A
        Vec u, v
        PetscErrorCode ierr

        !...local variable
        PetscScalar, Pointer, DIMENSION(:) :: uu, vv !...pointer to Petsc vector u/v, note: uu and vv are of size m_nOwnedDofs
        PetscErrorCode p_err
        INTEGER(int_p) :: i,b
        INTEGER(int_p) :: batchId, idxOffset_local, idxOffset_owned

        !...point to petsc vector v allowing vv to read and write to v
        CALL VecGetArrayF90(v, vv, p_err); CHKERRA(p_err)

        !...point to petsc vector u allowing uu to (only) read v
        CALL VecGetArrayReadF90(u, uu, p_err); CHKERRA(p_err)

        !...put local-size uu to ghost-included-size m_u (uu is 1-based index)
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                m_u(m_lcSid2LcMap(i,b)) = uu(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO

        !...matvec on ghost-included-size m_v = A * m_u
        !CALL profiler_start(m_mvTimer)
        CALL matvec_ghosted_bch(m_v, m_u)
        !CALL profiler_stop(m_mvTimer)

        !...copy m_v back to local-size array vv (vv is 1-based index)
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                vv(m_ocSid2OcMap(i,b)) = m_v(m_lcSid2LcMap(i,b))
            ENDDO
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(v, vv, p_err);
    END SUBROUTINE MatMult_amat_bch

    !...batched version of MatGetDiagonal_amat
    SUBROUTINE MatGetDiagonal_amat_bch(A, v, ierr)

        IMPLICIT NONE
        
        Mat A
        Vec v
        PetscErrorCode ierr

        !...local variable
        PetscScalar, Pointer, DIMENSION(:) :: vv !...pointer to Petsc vector v, note: vv of size m_nOwnedDofs
        PetscErrorCode p_err
        INTEGER(int_p) :: i,b
        INTEGER(int_p) :: batchId, idxOffset_local, idxOffset_owned

        !...point to petsc vector v allowing vv to read and write to v
        CALL VecGetArrayF90(v, vv, p_err); CHKERRA(p_err)

        !...set diagonal terms to PETSc vector
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                vv(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(v, vv, p_err);
    END SUBROUTINE MatGetDiagonal_amat_bch

    !...batched version of build (local owned) diagonal block petsc matrix to be used in bjacobi/sor preconditioner
    SUBROUTINE build_petsc_diagonal_block_matrix_bch
        IMPLICIT NONE

        !...local variable
        INTEGER(int_p) :: i, j, count, currf, b
        PetscScalar, DIMENSION(m_maxFacesPerCell+1) :: values
        PetscInt, DIMENSION(m_maxFacesPerCell+1) :: colInds
        PetscInt rowIdx
        
        !...allocate PETSc matrix for diagonal block
        IF (.not.(m_matPC_flag)) THEN
            CALL MatCreate(m_comm, mp_matPC, mp_ierr); CHKERRA(mp_ierr)
            CALL MATSetSizes(mp_matPC, m_nOwnedDofs_bch, m_nOwnedDofs_bch, PETSC_DECIDE, PETSC_DECIDE, mp_ierr); CHKERRA(mp_ierr)
            CALL MatSetFromOptions(mp_matPC, mp_ierr); CHKERRA(mp_ierr)
            IF (m_numprocs > 1) THEN
                CALL MatSetType(mp_matPC, MATMPIAIJ, mp_ierr); CHKERRA(mp_ierr)
                CALL MatMPIAIJSetPreallocation(mp_matPC, m_maxFacesPerCell + 1, PETSC_NULL_INTEGER, m_maxFacesPerCell, PETSC_NULL_INTEGER, mp_ierr); CHKERRA(mp_ierr)
            ELSE
                CALL MatSetType(mp_matPC, MATSEQAIJ, mp_ierr); CHKERRA(mp_ierr)
                CALL MatSeqAIJSetPreallocation(mp_matPC, m_maxFacesPerCell + 1, PETSC_NULL_INTEGER, mp_ierr); CHKERRA(mp_ierr)
            ENDIF
            m_matPC_flag = .TRUE.
        ENDIF

        !...set the diagonal block matrix
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                !...diagonal component is owned by my rank
                count = 1
                rowIdx = m_ocSid2GcMap(i,b) - 1 !...global index associated with batch b and owned cell i (no-batch version = lc_to_pgc(i) - 1)
                !...the first term is always the diagonal
                values(count) = m_matDiag(m_ocSid2OcMap(i,b))
                colInds(count) = rowIdx
                !...off-diagonal components
                DO j = 1,nfcell(i)
                    currf = lcf(i,j)
                    !...off-diagonal terms exist only when face lcf(i,j) is not a boundary face
                    IF (bface(currf) == 0) THEN
                        !...only set the components belong to my rank
                        IF ((lc_to_pgc(lcc(j,i)) .ge. m_globalDofBegin).and.(lc_to_pgc(lcc(j,i)) .le. m_globalDofEnd)) THEN
                            count = count + 1
                            values(count) = m_matOffDiag(j,m_ocSid2OcMap(i,b))
                            !...global index associated with batch b and OWNED cell lcc(j,i) because we are inside the if condition that LOCAL cell lcc(j,i) owned by me 
                            !...(no-batch version = lc_to_pgc(lcc(j,i)) - 1)
                            colInds(count) = m_ocSid2GcMap(lcc(j,i),b) - 1
                        ENDIF
                    ENDIF
                ENDDO
                CALL MatSetValues(mp_matPC, 1, rowIdx, count, colInds, values, INSERT_VALUES, mp_ierr)
            ENDDO
        ENDDO

        !...is there an alternative for sequential matrix?
        CALL MatAssemblyBegin(mp_matPC, MAT_FINAL_ASSEMBLY, mp_ierr)
        CALL MatAssemblyEnd(mp_matPC, MAT_FINAL_ASSEMBLY, mp_ierr)

    END SUBROUTINE build_petsc_diagonal_block_matrix_bch

    !...this is completely the same as MatGetDiagonalBlock_wrapper
    SUBROUTINE MatGetDiagonalBlock_wrapper_bch(A, aa, ierr)
        IMPLICIT NONE
        Mat A
        Mat aa
        PetscErrorCode ierr
        
        CALL MatGetDiagonalBlock(mp_matPC, aa, mp_ierr)

    END SUBROUTINE MatGetDiagonalBlock_wrapper_bch

    !...this is completely the same as MatSOR_wrapper
    SUBROUTINE MatSOR_wrapper_bch(A, b, omega, flag, shift, its, lits, x, ierr)
        IMPLICIT NONE
        
        Mat A
        Vec b
        PetscReal omega
        MatSORType flag
        PetscReal shift
        PetscInt its
        PetscInt lits
        Vec x
        PetscErrorCode ierr

        !CALL profiler_start(m_precondTimer)
        CALL MatSOR(mp_matPC, b, omega, flag, shift, its, lits, x, ierr)
        !CALL profiler_stop(m_precondTimer)
    END SUBROUTINE MatSOR_wrapper_bch

    !...batched version of matvec_ghosted
    SUBROUTINE matvec_ghosted_bch(v, u)
        USE GRID, ONLY : cell_type
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(inout) :: u
        REAL(real_p), DIMENSION(:), INTENT(out) :: v

        !...initialize v will be carried out in partial_matvec

        !...send values of u at owned cells to others who need as ghost, receive values of u at ghost cells from their owners
        CALL ghost_receive_begin_bch(u)

        !...while waiting for communication, perform matvec on independent cells
        CALL partial_matvec_bch(v, u, m_indCells, m_nIndCells)

        !...finishing communication
        CALL ghost_receive_end_bch(u)

        !...once communication is completed, perform matvec on dependent elements
        CALL partial_matvec_bch(v, u, m_depCells, m_nDepCells)
    END SUBROUTINE matvec_ghosted_bch


    !...batched version of partial_matvec
    SUBROUTINE partial_matvec_bch(v, u, cellList, listSize)
        USE GRID, ONLY : cell_type
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(in) :: u
        REAL(real_p), DIMENSION(:), INTENT(out) :: v
        INTEGER(int_p), DIMENSION(:), INTENT(in) :: cellList    !...list of cells to be performed v = K * u
        INTEGER(int_p), INTENT(in) :: listSize                  !...size of cellList

        !...local variables
        INTEGER(int_p) :: i, j, k, b

        DO b = 1,m_batchCount
            !...loop over cell list
            DO k = 1,listSize
                i = cellList(k) !...cell (local) id
                v(m_lcSid2LcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) * u(m_lcSid2LcMap(i,b)) !...diagonal term, this also initializes
                !...loop over off-diagonal terms
                DO j = 1,nfcell(i)
                    v(m_lcSid2LcMap(i,b)) = v(m_lcSid2LcMap(i,b)) + m_matOffDiag(j, m_ocSid2OcMap(i,b)) * u(m_lcSid2LcMap(lcc(j,i),b))
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE partial_matvec_bch


    !...batched version of set_aMat_matrix_n_rhs
    !...storing rhs in m_rhs first, and set to petsc vector mp_rhs at once
    SUBROUTINE compute_aMat_matrix_n_rhs_bch
        IMPLICIT NONE

        !...local variables
        !PetscErrorCode ierr   !...error code for PETSc subroutines
        INTEGER(int_p) :: i, j, mat, f, currf, ibface, bcon, pol, sout, mat1
        INTEGER(int_p) :: ic_1, ic_2, pol1, pol2
        INTEGER(int_p) :: b

        REAL(real_p) :: t_12, r_12, t_21, r_21, jface1, jface2, jface11, jface22
        REAL(real_p) :: velener_de, velener_r, n_omega1, n_omega2, disfun1

        REAL(real_p) :: temper, vol, gamma, sdotn, insdotn, vomega, vel, alpha, &
            rmuval, rxival, retval, inrmuval, inrxival, inretval, tempwall, &
            actin, actout, beta, mfp, jout, temper1

        !...initialize
        m_matDiag(:) = zero
        m_matOffDiag(:,:) = zero
        m_rhs(:) = zero
        !...loop over systems in the batch
        DO b = 1,m_batchCount
            !...get context data for system batch
            pol = polar_type(iband_ctx, p_ctx_bch(b))
            vel = gpvel(iband_ctx, pol)
            vomega = omega(si_ctx_bch(b))
            rmuval = rmu(si_ctx_bch(b))
            rxival = rxi(si_ctx_bch(b))
            retval = ret(si_ctx_bch(b))
            inrmuval = inrmu(si_ctx_bch(b))
            inrxival = inrxi(si_ctx_bch(b))
            inretval = inret(si_ctx_bch(b))

            !...loop over cells in the list
            DO i = 1, m_nOwnedDofs
                mat = cell_type(i)          !! material_type
                !IF (mat == ALUMINUM) CYCLE
                temper = tnot(i)

                !f = checkpolar(mat, iband_ctx, p_ctx_bch(b))
                !IF (f == 0) GO TO 100

                vol = volcell(i) * vomega

                CALL relaxtime(pol, mat, iband_ctx, temper, beta) !... what beta ?
                ! beta = 1.0D-6
                mfp = vel/beta

                !...global index for petsc vector rhs (even for CONJ case: ic_1 = i or ic_2 = i)
                m_idxs(m_ocSid2OcMap(i,b)) = m_ocSid2GcMap(i,b) - 1
                
                DO j = 1,nfcell(i)
                    currf = lcf(i,j) !... face id of face j of cell i
                    gamma = mfp * areaf(currf)/vol
                    
                    IF (threed) THEN
                        !... sdotn is the dot product of si_ctx (rmuval) and nf (vecf), what is insdotn?
                        sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival + vecfz(currf)*retval        
                        insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival+ vecfz(currf)*inretval
                    ELSE 
                        sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival     
                        insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival
                    ENDIF
                    ! sdotn = sdotn_array(currf, si_ctx_bch(b))
                    ! insdotn = insdotn_array(currf, si_ctx_bch(b))
                    
                    IF (bface(currf) == 0) THEN
                        !...Interior face, no contribution to rhs
                        alpha = normdir(i,j)
                        sdotn = alpha * sdotn
                        insdotn = alpha*insdotn

                        m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))
                        m_matOffDiag(j, m_ocSid2OcMap(i,b)) = -Max(0.0,-sdotn)*gamma*((insdotn/(sdotn+tiny)))

                    ELSE
                        !... Boundary face
                        ibface = f_to_bf(currf)
                        bcon = bctype(ibface)
                        tempwall = temp_bc(ibface)
                        
                        CALL iwalls(pol, mat, iband_ctx, tempwall, actin)
                        
                        jout = zero 
                        
                        SELECT CASE(bcon)
                        CASE(ADIA)
                            jout = jfacepos(ibface, pol, iband_ctx)

                            !... find the reflection direction
                            CALL reflec(si_ctx_bch(b), currf, sdotn, sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)

                            !...store matrix
                            m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                            !...store rhs vector
                            m_rhs(m_ocSid2OcMap(i,b)) = m_rhs(m_ocSid2OcMap(i,b)) + &
                                    dos*intensity(i, sout, icount_ctx, pol) * gamma * &
                                    Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) &
                                    +((one-dos)*jout) * gamma * Max(0.0,-sdotn)* &
                                    ((insdotn/(sdotn+tiny)))/pi

                        CASE(ISOT)
                            !...store matrix
                            m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                            !... store rhs vector
                            !rowIdx = m_ocSid2GcMap(i,b) - 1
                            m_rhs(m_ocSid2OcMap(i,b)) = m_rhs(m_ocSid2OcMap(i,b)) + &
                                actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                        CASE (HFLUX) !!! assuming a virtual block of aluminium!!
                            !...store matrix
                            m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                            !... store rhs vector
                            m_rhs(m_ocSid2OcMap(i,b)) = m_rhs(m_ocSid2OcMap(i,b)) + &
                                actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                        CASE (ALSI)
                            !...store matrix
                            m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                            !... store rhs vector
                            m_rhs(m_ocSid2OcMap(i,b)) = m_rhs(m_ocSid2OcMap(i,b)) + &
                                actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))
                        
                        !!! symmetry boundary
                        CASE(SYMM)
                            !...store matrix
                            jout = jfacepos(ibface, pol, iband_ctx)
                            CALL reflec(si_ctx_bch(b), currf, sdotn, sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)
                            
                            m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                            !... store rhs vector
                            m_rhs(m_ocSid2OcMap(i,b)) = m_rhs(m_ocSid2OcMap(i,b)) + &
                                intensity(i,sout,icount_ctx,pol) * gamma * &
                                            Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))
                           
                        !!! INTERFACE....     
                        ! CASE (CONJ)
                        !     ic_1 = lfc(currf, 1)
                        !     ic_2 = lfc(currf, 2)
                            
                        !     pol1 = pol
                        !     pol2 = pol
                            
                        !     IF (cell_type(ic_1) /= cell_type(ic_2)) THEN
                        !         pol1 = mat_polar(cell_type(ic_1),pol)
                        !         pol2 = mat_polar(cell_type(ic_2),pol)
                        !     ENDIF
                        
                        !     CALL trans_coeff (iband_ctx, pol1, pol2, ic_1, ic_2, ibface, t_12)

                        !     !!! file name changed, compilation done according to need              
                            
                        !     !t_12 = 0.5
                        !     r_12 = 1 - t_12
                        !     t_21 = 1 - t_12
                        !     !t_21 = t_12
                        !     r_21 = 1 - t_21
                                                    
                        !     jface1 = jfacepos(ibface, pol, iband_ctx)
                        !     jface2 = jfaceneg(ibface, pol, iband_ctx)
                        !     jface11 = jfacepos(ibface, pol1, iband_ctx)
                        !     jface22 = jfaceneg(ibface, pol2, iband_ctx)

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
                        !         !...store matrix
                        !         m_matDiag(m_ocSid2OcMap(ic_1,b)) = m_matDiag(m_ocSid2OcMap(ic_1,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !         !... store rhs vector
                        !         !rowIdx = m_ocSid2GcMap(ic_1,b) - 1
                        !         m_rhs(m_ocSid2OcMap(ic_1,b)) = m_rhs(m_ocSid2OcMap(ic_1,b)) + &
                        !             velener_r*((jface1*r_12)/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) + &
                        !                     (jface22*t_21/pi)*gamma*Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))*velener_r

                        !     ELSE
                        !         !...store matrix
                        !         m_matDiag(m_ocSid2OcMap(ic_2,b)) = m_matDiag(m_ocSid2OcMap(ic_2,b)) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !         !... store rhs vector
                        !         !rowIdx = m_ocSid2GcMap(ic_2,b) - 1
                        !         m_rhs(m_ocSid2OcMap(ic_2,b)) = m_rhs(m_ocSid2OcMap(ic_2,b)) + &
                        !             velener_r*((jface2*r_21)/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) + &
                        !                     (jface11*t_12/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))*velener_r

                        !     ENDIF                                                                  
                        END SELECT !! case of bcon (type of BCs of face j of cell i)
                    ENDIF !Bounadry or Interior Face

                    !WRITE (*,*)"bsnord 4", iband_ctx, pol, si_ctx, i,j
                    !  scs_tran(i,pol) = scs_tran(i,pol)*gamma       !!! added 28th aug 
                ENDDO ! Faces of cell    loop j
                    !WRITE (*,*)"bsnord 5", iband_ctx, pol, si_ctx, i
                    !!!!!13thaug$$$!!!!!!!!!
            100 ENDDO       !! 1st cell loop 

            !... now compute the contribution of the first and the third terms of eq (3.17)
            DO i = 1,m_nOwnedDofs
                mat1 = cell_type(i)
                ! IF  (mat1== ALUMINUM) CYCLE 
                temper1= tnot(i)
                !scs(i) = scs(i) + scs_tran(i,pol)       !!! added 28th aug 
                !f = checkpolar(mat1, iband_ctx, p_ctx_bch(b))
                !IF (f /= 0) THEN
                    CALL relaxtime(pol, mat1, iband_ctx, temper1, beta)
                    !beta = 1.0D-6
                    CALL iwalls(pol, mat1, iband_ctx, temper1, actout)

                    !...store matrix
                    m_matDiag(m_ocSid2OcMap(i,b)) = m_matDiag(m_ocSid2OcMap(i,b)) + (idt/beta + one)

                    !...store rhs vector
                    m_rhs(m_ocSid2OcMap(i,b)) = m_rhs(m_ocSid2OcMap(i,b)) + &
                        actout + idt *intensityone(i, si_ctx_bch(b), icount_ctx, pol)/beta
                !ENDIF
            ENDDO
        ENDDO
    END SUBROUTINE compute_aMat_matrix_n_rhs_bch

    !...this version includes computing the rhs norms
    SUBROUTINE modify_aMat_rhs_for_dx_solving_bch_1
        IMPLICIT NONE

        !...local variables
        INTEGER(int_p) :: i, b, si, p
        REAL(real_p) :: rhsValue

        !...put intensity of previous iteration into m_u (we used m_u for this subroutine to avoid allocation of new memory)
        DO b = 1,m_batchCount
            si = si_ctx_bch(b)
            p = p_ctx_bch(b)
            DO i = 1,m_nOwnedDofs
                m_u(m_lcSid2LcMap(i,b)) = intensity(i, si, icount_ctx, p)
            ENDDO
        ENDDO

        !...peform m_v = A * m_u where A is the matrix of the system of equations of intensity of current iteration
        CALL matvec_ghosted_bch(m_v, m_u)

        !...subtract m_v from the m_rhs
        DO b = 1,m_batchCount
            m_rhsNorm_bch(b) = 0.d0
            DO i = 1,m_nOwnedDofs
                rhsValue = m_rhs(m_ocSid2OcMap(i,b)) - m_v(m_lcSid2LcMap(i,b))
                m_rhs(m_ocSid2OcMap(i,b)) = rhsValue;
                m_rhsNorm_bch(b) = m_rhsNorm_bch(b) + (rhsValue * rhsValue)
            ENDDO
        ENDDO
        CALL MPI_ALLREDUCE(m_rhsNorm_bch(1), m_rhsNormSum_bch(1), m_batchCount, MPI_REAL8, MPI_SUM, m_comm, mp_ierr)
        DO b = 1,m_batchCount
            m_rhsNormSum_bch(b) = SQRT(m_rhsNormSum_bch(b))
        ENDDO
    END SUBROUTINE modify_aMat_rhs_for_dx_solving_bch_1


    SUBROUTINE compute_residual_bch(sol, res)
        IMPLICIT NONE
        Vec, INTENT(IN) :: sol
        REAL(real_p), DIMENSION(:), INTENT(OUT) :: res !...array of residuals of each linear system in batch

        !...local variables
        PetscScalar, Pointer, DIMENSION(:) :: ssol
        PetscErrorCode p_err
        INTEGER(int_p) :: i,b

        REAL(real_p), DIMENSION(:), ALLOCATABLE :: temp1, temp2, temp1_sum, temp2_sum
        INTEGER(int_p), PARAMETER :: root_proc = 0

        ALLOCATE(temp1(m_batchCount), temp2(m_batchCount))
        ALLOCATE(temp1_sum(m_batchCount), temp2_sum(m_batchCount))

        !... transfer (petsc) vector to ghosted m_u
        CALL VecGetArrayReadF90(sol, ssol, p_err)
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                m_u(m_lcSid2LcMap(i,b)) = ssol(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO

        !...perform m_v = A * m_u
        CALL matvec_ghosted_bch(m_v, m_u)

        !...subtract rhs from m_v
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                m_v(m_lcSid2LcMap(i,b)) = m_v(m_lcSid2LcMap(i,b)) - m_rhs(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO

        !...compute L2 norm of m_v and rhs
        DO b = 1,m_batchCount
            temp1(b) = 0.d0 
            temp2(b) = 0.d0
            DO i = 1,m_nOwnedDofs
                temp1(b) = temp1(b) + m_v(m_lcSid2LcMap(i,b)) * m_v(m_lcSid2LcMap(i,b))
                temp2(b) = temp2(b) + m_rhs(m_ocSid2OcMap(i,b)) * m_rhs(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO
        
        !...reduction sum to root
        DO b = 1,m_batchCount
            CALL MPI_REDUCE(temp1(b), temp1_sum(b), 1, MPI_REAL8, MPI_SUM, root_proc, m_comm, p_err)
            CALL MPI_REDUCE(temp2(b), temp2_sum(b), 1, MPI_REAL8, MPI_SUM, root_proc, m_comm, p_err)
        ENDDO

        !...root computes norm of residual
        IF (m_rank .eq. root_proc) THEN
            DO b = 1,m_batchCount
                temp1_sum(b) = SQRT(temp1_sum(b))
                temp2_sum(b) = SQRT(temp2_sum(b))
                res(b) = temp1_sum(b)/temp2_sum(b)
            ENDDO
        ENDIF

        !...broadcast res computed by root to other ranks
        DO b = 1,m_batchCount
            CALL MPI_BCAST(res(b), 1, MPI_REAL8, root_proc, m_comm, p_err)
        ENDDO

        DEALLOCATE(temp1, temp1_sum, temp2, temp2_sum)

    END SUBROUTINE compute_residual_bch


    !...set aMat rhs to PETSc vector
    SUBROUTINE set_petsc_rhs_bch(rhs)
        IMPLICIT NONE
        Vec rhs

        CALL VecSetValues(rhs, m_nOwnedDofs_bch, m_idxs, m_rhs, INSERT_VALUES, mp_ierr)
        !...assemble petsc rhs vector
        CALL VecAssemblyBegin(mp_rhs, mp_ierr)
        CALL VecAssemblyEnd(mp_rhs, mp_ierr)
    END SUBROUTINE set_petsc_rhs_bch


    !...set s to the intensity of previous iteration, to be used as non-zero initial guess
    !...batched version
    SUBROUTINE set_sol_prev_itsy_bch(s)
        IMPLICIT NONE
        Vec s

        !...internal variables
        INTEGER(int_p) :: i, j, si, p
        PetscScalar, Pointer, DIMENSION(:) :: ss
        PetscErrorCode p_err

        !...point to petsc vector allowing ss to read/write to sol
        CALL VecGetArrayF90(s, ss, p_err); CHKERRA(p_err)

        !...set solution of batched system (used as non-zero initial guess) from intensity array
        DO j = 1,m_batchCount
            si = si_ctx_bch(j)
            p = p_ctx_bch(j)
            DO i = 1,m_nOwnedDofs
                ss(m_ocSid2OcMap(i,j)) = intensity(i, si, icount_ctx, p)
            ENDDO
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(s, ss, p_err); CHKERRA(p_err)
        
    END SUBROUTINE set_sol_prev_itsy_bch


    !...compute PETSc global rhs vector (only used to compute norm in batch method)
    SUBROUTINE compute_norm_rhs_dir(norm)
        IMPLICIT NONE
        REAL(real_p), INTENT(OUT) :: norm

        !...local variables
        !PetscErrorCode ierr   !...error code for PETSc subroutines
        INTEGER(int_p) :: i, j, mat, currf, ibface, bcon, pol, sout, mat1
        INTEGER(int_p) :: ic_1, ic_2, pol1, pol2

        REAL(real_p) :: temper, vol, gamma, sdotn, insdotn, vomega, vel, alpha, &
            rmuval, rxival, retval, inrmuval, inrxival, inretval, tempwall, &
            actin, actout, beta, mfp, jout, temper1

        PetscInt rowIdx
        PetscScalar :: rhs_value

        pol = polar_type(iband_ctx, p_ctx)
        vel = gpvel(iband_ctx, pol)
        
        vomega = omega(si_ctx)
        rmuval = rmu(si_ctx)
        rxival = rxi(si_ctx)
        retval = ret(si_ctx)
        inrmuval = inrmu(si_ctx)
        inrxival = inrxi(si_ctx)
        inretval = inret(si_ctx)

        !...initialize norm
        norm = 0.d0

        !... main loop finds the coefficients the rhs b
        DO i = 1, m_nOwnedDofs
            mat = cell_type(i)          !! material_type
            !IF (mat == ALUMINUM) CYCLE
            temper = tnot(i)

            !f = checkpolar(mat,iband_ctx, p_ctx)
            !IF (f == 0) GO TO 100

            vol = volcell(i) * vomega

            CALL relaxtime(pol, mat, iband_ctx, temper, beta) !... what beta ?
            ! beta = 1.0D-6
            mfp = vel/beta
            !write(*,'(a,i0,a,f,a,f,a,i0)')'i=',i,', vel=',vel,', beta=',beta
            
            rowIdx = lc_to_pgc(i) - 1 !... row index = global id of cell i (note PETSc uses 0-based index for both C and Fortran)
            
            rhs_value = 0.0 !...value of rhs corresponding to cell i

            DO j = 1,nfcell(i)
                currf = lcf(i,j) !... face id of face j of cell i
                gamma = mfp * areaf(currf)/vol

                IF (threed) THEN
                    !... sdotn is the dot product of Si (rmuval) and nf (vecf), what is insdotn?
                    sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival + vecfz(currf)*retval        
                    insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival+ vecfz(currf)*inretval
                ELSE 
                    sdotn = vecfx(currf) * rmuval + vecfy(currf) * rxival     
                    insdotn = vecfx(currf) * inrmuval + vecfy(currf) * inrxival
                ENDIF
                ! sdotn = sdotn_array(currf, si_ctx)
                ! insdotn = insdotn_array(currf, si_ctx)
                
                IF (bface(currf) .eq. 1) THEN
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

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + dos*intensity(i,sout,icount_ctx,pol) * gamma * &
                                Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) &
                                +((one-dos)*jout) * gamma * Max(0.0,-sdotn)* &
                                ((insdotn/(sdotn+tiny)))/pi

                    CASE(ISOT)
                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    CASE (HFLUX) !!! assuming a virtual block of aluminium!!
                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    CASE (ALSI)
                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))
                    
                    !!! symmetry boundary
                    CASE(SYMM)
                        jout = jfacepos(ibface,pol,iband_ctx)
                        CALL reflec(si_ctx,currf,sdotn,sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)

                        !... assemble PETSc rhs vector
                        rhs_value = rhs_value + intensity(i,sout,icount_ctx,pol) * gamma * &
                                        Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    END SELECT !! case of bcon (type of BCs of face j of cell i)
                ENDIF !Bounadry or Interior Face  
                !  scs_tran(i,pol) = scs_tran(i,pol)*gamma       !!! added 28th aug 
            ENDDO ! Faces of cell    loop j
          
            !...accumulate rhs of row/cell i to norm
            norm = norm + (rhs_value * rhs_value)

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
                CALL iwalls(pol, mat1, iband_ctx, temper1, actout)
                rhs_value = actout + idt *intensityone(i,si_ctx,icount_ctx,pol)/beta
                norm = norm + (rhs_value * rhs_value)
            !ENDIF
        ENDDO

        !...take square-root of sum
        norm = SQRT(norm)

    END SUBROUTINE compute_norm_rhs_dir


    !...user-defined convergence test
    SUBROUTINE MyKSPConverged(ksp, n, rnorm, flag, dummy, ierr)
        IMPLICIT NONE
        KSP ksp
        PetscErrorCode ierr
        PetscInt n, dummy
        KSPConvergedReason flag
        PetscReal rnorm

        !...internal variables
        Vec sol_bch
        PetscMPIInt myrank
        INTEGER(int_p) :: i,b
        PetscScalar, Pointer, DIMENSION(:) :: sol_v

        !...get intermediate solution
        CALL KSPBuildSolution(ksp, PETSC_NULL_VEC, sol_bch, mp_ierr)

        !...point to intermediate solution
        CALL VecGetArrayF90(sol_bch, sol_v, mp_ierr); CHKERRA(mp_ierr)

        !...transfer petsc vector to (ghosted) m_ui (here I use m_ui and m_vi to not "touch" the m_u and m_v in matvec)
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                m_ui(m_lcSid2LcMap(i,b)) = sol_v(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO

        !...perform m_vi = A * m_ui
        CALL matvec_ghosted_bch(m_vi, m_ui)

        !...subtract rhs from m_vi to get the residual vector of the batched system
        DO b = 1,m_batchCount
            DO i = 1,m_nOwnedDofs
                !...this is residual vector (portion that I own) of batched system
                m_vi(m_lcSid2LcMap(i,b)) = m_vi(m_lcSid2LcMap(i,b)) - m_rhs(m_ocSid2OcMap(i,b))
            ENDDO
        ENDDO

        !...compute norms of each system in the batch and divide by the norm of rhs of each system
        DO b = 1,m_batchCount
            !...compute sum of square of residual_i (portion that I own)
            m_resNorm_bch(b) = 0.d0
            DO i = 1,m_nOwnedDofs
                m_resNorm_bch(b) = m_resNorm_bch(b) + (m_vi(m_lcSid2LcMap(i,b)) * m_vi(m_lcSid2LcMap(i,b)))
            ENDDO
        ENDDO
        !...communicate with others to have the total sum across ranks
        CALL MPI_ALLREDUCE(m_resNorm_bch(1), m_resNormSum_bch(1), m_batchCount, MPI_REAL8, MPI_SUM, m_comm, mp_ierr)
        !...take square root of sum of square
        DO b = 1,m_batchCount
            m_resNormSum_bch(b) = SQRT(m_resNormSum_bch(b)) !...this is norm of absolute residual
        ENDDO

        !...check the convergence:
        flag = 1
        DO b = 1,m_batchCount
            IF (m_resNormSum_bch(b) > MAX(m_rhsNormSum_bch(b) * mp_RTOL, mp_ATOL)) THEN
                flag = 0
                EXIT
            ENDIF
        ENDDO
        
    END SUBROUTINE MyKSPConverged
END MODULE amat_batched