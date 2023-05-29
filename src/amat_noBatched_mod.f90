! \brief module to solve using assembly-free aMat no batching
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! Solve the linear systems of intensity using assembly-free method, no batching
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
MODULE amat_noBatched
    USE petsc_solver
    USE amat_parent

    IMPLICIT NONE

    CONTAINS
    
    !...initialize aMat variables
    SUBROUTINE init_aMat
        IMPLICIT NONE

        INTEGER(int_p) :: ierr

        !...local-size (including ghost cells) vectors u and v used in matvec function
        ALLOCATE(m_u(m_nLocalDofs), stat=ierr)
        ALLOCATE(m_v(m_nLocalDofs), stat=ierr)

        !...diagonal terms of rows that I own
        ALLOCATE(m_matDiag(m_nOwnedDofs), stat=ierr)

        !...determine max number of faces per cell
        CALL determine_max_faces_per_cell(m_maxFacesPerCell)

        !...off-diagonal terms of rows that I own
        ALLOCATE(m_matOffDiag(m_maxFacesPerCell, m_nOwnedDofs), stat=ierr) !...Fortran uses column-major order (consecutive elements of a column are stored next to each other)

        !...rhs vector, used to store rhs during assemble and only set to Petsc vector once (added on Oct 18 2021)
        ALLOCATE(m_rhs(m_nOwnedDofs), stat=ierr)
        ALLOCATE(m_idxs(m_nOwnedDofs), stat=ierr) !...array of global cell ids

        !...determine dep/ind cells
        CALL identifyIndCells

        !...profiling
        !CALL profiler_init(m_mvTimer)
        !CALL profiler_init(m_precondTimer)
        
        !...flag if block-jacobi matrix is allocated
        m_matPC_flag = .FALSE.

    END SUBROUTINE init_aMat
    
    !...subroutine to finalize aMat
    SUBROUTINE finalize_aMat
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
        IF (ALLOCATED(m_statusArray)) DEALLOCATE(m_statusArray)
        IF (ALLOCATED(m_u)) DEALLOCATE(m_u)
        IF (ALLOCATED(m_v)) DEALLOCATE(m_v)
        IF (ALLOCATED(m_indCells)) DEALLOCATE(m_indCells)
        IF (ALLOCATED(m_depCells)) DEALLOCATE(m_depCells)
        IF (ALLOCATED(m_matDiag)) DEALLOCATE(m_matDiag)
        IF (ALLOCATED(m_matOffDiag)) DEALLOCATE(m_matOffDiag)
        IF (ALLOCATED(m_rhs)) DEALLOCATE(m_rhs) !...Oct 18 2021
        IF (ALLOCATED(m_idxs)) DEALLOCATE(m_idxs) !...Oct 18 2021
        IF (m_matPC_flag) THEN
            CALL MatDestroy(mp_matPC, mp_ierr); CHKERRA(mp_ierr)
        ENDIF
        CALL MatDestroy(mp_mat, mp_ierr) !...although we create mp_mat as a shell, we should destroy it
        IF (ALLOCATED(m_sorX)) DEALLOCATE(m_sorX)

    END SUBROUTINE finalize_aMat

    !..."allocate" matrix shell
    SUBROUTINE allocate_petsc_matshell
        IMPLICIT NONE

        !...create matrix shell
        CALL MatCreateShell(m_comm, m_nOwnedDofs, m_nOwnedDofs, PETSC_DECIDE, PETSC_DECIDE, PETSC_NULL_INTEGER, mp_mat, mp_ierr)

        !...define matrix-vector multiplication operator
        CALL MatShellSetOperation(mp_mat, MATOP_MULT, MatMult_amat, mp_ierr)

        !...define get matrix diagonal operator, used in pc jacobi
        CALL MatShellSetOperation(mp_mat, MATOP_GET_DIAGONAL, MatGetDiagonal_amat, mp_ierr)

        !...define get diagonal block matrix operator, used in pc bjacobi
        CALL MatShellSetOperation(mp_mat, MATOP_GET_DIAGONAL_BLOCK, MatGetDiagonalBlock_wrapper, mp_ierr)
        !CALL MatShellSetOperation(mp_mat, MATOP_GET_DIAGONAL_BLOCK, MatGetDiagonalBlock_amat, mp_ierr)

        !...define get sor sweep, used in pc sor
        CALL MatShellSetOperation(mp_mat, MATOP_SOR, MatSOR_wrapper, mp_ierr)
        !CALL MatShellSetOperation(mp_mat, MATOP_SOR, MatSOR_amat, mp_ierr)

    END SUBROUTINE allocate_petsc_matshell

    
    !...subroutine to build scatter map
    SUBROUTINE buildScatterMap
        IMPLICIT NONE

        !...local variables
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
        
        !...list of ghost ids in terms of petsc global id
        ALLOCATE(ghostGIds(m_nGhostDofs), stat=ierr)
        DO i = 1,m_nGhostDofs
            ghostGIds(i) = lc_to_pgc(m_nOwnedDofs + i)
            !write(*,*) 'rank ',m_rank,', ghostGId(',i,')=',ghostGIds(i)
        ENDDO

        !...determine the owner of ghost cells (required: global cell ids are continuously partitioned across processors)
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
        m_totalSend = nRecvs

        !...total cells to be received from
        m_totalRecv = m_nGhostDofs

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
            m_sendCellCounts(i) = recvCounts(i)
            m_sendCellOffset(i) = recvOffset(i)
            m_recvCellCounts(i) = ghostCounts(i)
            m_recvCellOffset(i) = ghostOffset(i)
        ENDDO
        DO i = 1,nRecvs
            globalId = recvBuff(i) !...global ID of cell that I need to send to other ranks
            localId = globalId - localCellOffset(m_rank) !...convert to local cell id which is 1-based index
            m_sendCellIds(i) = localId
        ENDDO

        !...allocate buffers for send and recv in matvec
        IF (m_totalSend > 0) THEN
            ALLOCATE(m_sendBuf(0:m_totalSend-1), stat=ierr)
        ENDIF
        IF (m_totalRecv > 0) THEN
            ALLOCATE(m_recvBuf(0:m_totalRecv-1), stat=ierr)
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
        ALLOCATE(m_statusArray(MPI_STATUS_SIZE, m_nProcsSend + m_nProcsRecv), stat=ierr)

        IF (ALLOCATED(localCellCounts)) DEALLOCATE(localCellCounts)
        IF (ALLOCATED(localCellOffset)) DEALLOCATE(localCellOffset)
        IF (ALLOCATED(ghostCounts)) DEALLOCATE(ghostCounts)
        IF (ALLOCATED(ghostOffset)) DEALLOCATE(ghostOffset)
        IF (ALLOCATED(recvCounts)) DEALLOCATE(recvCounts)
        IF (ALLOCATED(recvOffset)) DEALLOCATE(recvOffset)
        IF (ALLOCATED(ghostGIds)) DEALLOCATE(ghostGIds)
        IF (ALLOCATED(ghostOwner)) DEALLOCATE(ghostOwner)
        IF (ALLOCATED(recvBuff)) DEALLOCATE(recvBuff)
    END SUBROUTINE buildScatterMap
    

    !...subroutine performs sending cells that I own to other processes who need them as ghost
    !...and receiving ghost values from their owners
    SUBROUTINE ghost_receive_begin(vec)
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(IN) :: vec

        !...local variables
        INTEGER(int_p) :: i, r, counter
        INTEGER :: mpierr
        INTEGER(int_p), PARAMETER :: commTag = 1

        IF (m_numprocs == 1) RETURN

        !...send values of v at (owned) cells to other ranks who need them as ghost
        counter = 0
        IF (m_totalSend > 0) THEN
            !...put values of vec at positions indcated by m_sendCellIds to sending buffer
            DO i = 1,m_totalSend 
                m_sendBuf(i-1) = vec(m_sendCellIds(i)) !...m_sendBuf is 0-based index
            ENDDO
            
            !...send data to other processes
            DO i = 1,m_nProcsSend 
                r = m_sendRankIds(i) !...rank that I am going to send to
                counter = counter + 1
                CALL MPI_Isend(m_sendBuf(m_sendCellOffset(r)), m_sendCellCounts(r), MPI_REAL8, r, &
                                commTag, m_comm, m_request(counter), mpierr)
            ENDDO
        ENDIF

        !...receive values of v at ghost cells from their owners
        IF (m_totalRecv > 0) THEN
            DO i = 1,m_nProcsRecv
                r = m_recvRankIds(i) !...rank that I received from
                counter = counter + 1
                CALL MPI_Irecv(m_recvBuf(m_recvCellOffset(r)), m_recvCellCounts(r), MPI_REAL8, r, &
                                commTag, m_comm, m_request(counter), mpierr)
            ENDDO
        ENDIF
    END SUBROUTINE ghost_receive_begin

    !...subroutine to synchronize the isend and irecv started in ghost_receive_begin
    SUBROUTINE ghost_receive_end(vec)
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(OUT) :: vec

        !...local variables
        INTEGER(int_p) :: i
        INTEGER :: mpierr

        IF (m_numprocs == 1) RETURN

        !...synchronize (wait) all requests of isend and irecv have started
        CALL MPI_Waitall(m_nProcsRecv+m_nProcsSend, m_request, m_statusArray, mpierr)
        
        !...put the received values to the ghost position of vec
        DO i = 1,m_totalRecv
            vec(i+m_nOwnedDofs) = m_recvBuf(i-1) !...m_recvBuf is 0-based index
        ENDDO
    END SUBROUTINE ghost_receive_end

    !...this subroutine is the interface of PETSc function MatMult(A, u, v, ierr) that performs v = A * u
    !...to be used in MatShell set operation
    !...u and v are Petsc vectors that only hold the local dofs (i.e. owned dofs)
    SUBROUTINE MatMult_amat(A, u, v, ierr)
        IMPLICIT NONE
        
        Mat A
        Vec u, v
        PetscErrorCode ierr

        !...local variable
        PetscScalar, Pointer, DIMENSION(:) :: uu, vv !...pointer to Petsc vector u/v, note: uu and vv are of size m_nOwnedDofs
        PetscErrorCode p_err
        INTEGER(int_p) :: i

        !...point to petsc vector v allowing vv to read and write to v
        CALL VecGetArrayF90(v, vv, p_err); CHKERRA(p_err)

        !...point to petsc vector u allowing uu to (only) read v
        CALL VecGetArrayReadF90(u, uu, p_err); CHKERRA(p_err)

        !...put local-size uu to ghost-included-size m_u (uu is 1-based index)
        DO i = 1,m_nOwnedDofs
            m_u(i) = uu(i)
        ENDDO

        !...matvec on ghost-included-size m_v = A * m_u
        !CALL profiler_start(m_mvTimer)
        CALL matvec_ghosted(m_v, m_u)
        !CALL profiler_stop(m_mvTimer)

        !...copy m_v back to local-size array vv (vv is 1-based index)
        DO i = 1,m_nOwnedDofs
            vv(i) = m_v(i)
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(v, vv, p_err);
    END SUBROUTINE MatMult_amat


    !...this subroutine is the interface of PETSc function MatGetDiagonal(A, v, ierr)
    !...to be used in MatShell set operation
    !...v is Petsc vector that only holds the local dofs (i.e. owned dofs)
    SUBROUTINE MatGetDiagonal_amat(A, v, ierr)
        IMPLICIT NONE
        
        Mat A
        Vec v
        PetscErrorCode ierr

        !...local variable
        PetscScalar, Pointer, DIMENSION(:) :: vv !...pointer to Petsc vector v, note: size of vv= m_nOwnedDofs
        PetscErrorCode p_err
        INTEGER(int_p) :: i

        !...point to petsc vector v allowing vv to read and write to v
        CALL VecGetArrayF90(v, vv, p_err); CHKERRA(p_err)

        !...set diagonal terms to PETSc vector
        DO i = 1,m_nOwnedDofs
            vv(i) = m_matDiag(i)
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(v, vv, p_err);
    END SUBROUTINE MatGetDiagonal_amat


    !...build (local owned) diagonal block petsc matrix to be used in bjacobi/sor preconditioner
    SUBROUTINE build_petsc_diagonal_block_matrix
        IMPLICIT NONE

        !...local variable
        INTEGER(int_p) :: i, j, count, currf
        PetscScalar, DIMENSION(m_maxFacesPerCell+1) :: values
        PetscInt, DIMENSION(m_maxFacesPerCell+1) :: colInds
        PetscInt rowIdx
        
        !...allocate PETSc matrix for diagonal block
        IF (.not.(m_matPC_flag)) THEN
            CALL MatCreate(m_comm, mp_matPC, mp_ierr); CHKERRA(mp_ierr)
            CALL MATSetSizes(mp_matPC, m_nOwnedDofs, m_nOwnedDofs, PETSC_DECIDE, PETSC_DECIDE, mp_ierr); CHKERRA(mp_ierr)
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
        DO i = 1,m_nOwnedDofs
            !...diagonal component is owned by my rank
            count = 1
            rowIdx = lc_to_pgc(i) - 1
            values(count) = m_matDiag(i)
            colInds(count) = rowIdx
            !...off-diagonal components
            DO j = 1,nfcell(i)
                currf = lcf(i,j)
                !...off-diagonal terms exist only when face lcf(i,j) is not a boundary face
                IF (bface(currf) == 0) THEN
                    !...only set the components belong to my rank
                    IF ((lc_to_pgc(lcc(j,i)) .ge. m_globalDofBegin).and.(lc_to_pgc(lcc(j,i)) .le. m_globalDofEnd)) THEN
                        count = count + 1
                        values(count) = m_matOffDiag(j,i)
                        colInds(count) = lc_to_pgc(lcc(j,i)) - 1
                    ENDIF
                ENDIF
            ENDDO
            CALL MatSetValues(mp_matPC, 1, rowIdx, count, colInds, values, INSERT_VALUES, mp_ierr)
        ENDDO

        !...is there an alternative for sequential matrix?
        CALL MatAssemblyBegin(mp_matPC, MAT_FINAL_ASSEMBLY, mp_ierr)
        CALL MatAssemblyEnd(mp_matPC, MAT_FINAL_ASSEMBLY, mp_ierr)

    END SUBROUTINE build_petsc_diagonal_block_matrix


    !...this subroutine is the interface of PETSc function MatGetDiagonalBlock(A, a, ierr), to be used in MatShell set operation
    SUBROUTINE MatGetDiagonalBlock_wrapper(A, aa, ierr)
        IMPLICIT NONE
        
        Mat A
        Mat aa
        PetscErrorCode ierr
        
        CALL MatGetDiagonalBlock(mp_matPC, aa, mp_ierr)
        
    END SUBROUTINE MatGetDiagonalBlock_wrapper


    !...this subroutine is the interface of PETSc function MatGetDiagonalBlock(A, a, ierr), to be used in MatShell set operation
    SUBROUTINE MatSOR_wrapper(A, b, omega, flag, shift, its, lits, x, ierr)
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
    END SUBROUTINE MatSOR_wrapper


    SUBROUTINE MatSOR_amat(A, b, omega, flag, shift, its, lits, x, ierr)
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

        !...local variable
        INTEGER(int_p) :: i,j,k,currf, err
        PetscScalar, Pointer, DIMENSION(:) :: xx !...pointer to Petsc vector x (note: size of xx = m_nOwnedDofs)
        PetscScalar, Pointer, DIMENSION(:) :: bb !...pointer to Petsc vector b (note: size of bb = m_nOwnedDofs)
        PetscErrorCode p_err

        PetscInt rowIdx
        PetscInt colIdx

        !...allocate vector for storing current-step x
        IF (.not.(ALLOCATED(m_sorX))) THEN
            AlLOCATE(m_sorX(m_nOwnedDofs), stat=err)
        ENDIF

        !...point to PETSc vector x, allowing xx to read and write to x
        CALL VecGetArrayF90(x, xx, p_err); CHKERRA(p_err)
        CALL VecGetArrayF90(b, bb, p_err); CHKERRA(p_err)

        !...loop over local cells
        DO i = 1,m_nOwnedDofs
            rowIdx = lc_to_pgc(i) !...no need to convert to 0-based index because we just compare relatively to column index
            
            !...rhs term and diagonal term
            m_sorX(i) = (omega/m_matDiag(i)) * bb(i) + (1.0 - omega) * xx(i)
            
            !...loop over all faces of cell i
            DO j = 1, nfcell(i)
                currf = lcf(i,j)
                !...off-diagonal term exists only when face lcf(i,j) is not a boundary face
                IF (bface(currf) == 0) THEN
                    colIdx = lc_to_pgc(lcc(j,i)) !...global index of column associated with face j
                    !...similar to block jacobi, only consider the terms belong to my rank
                    IF ((colIdx .ge. m_globalDofBegin).and.(colIdx .le. m_globalDofEnd)) THEN
                        !...upper part --> use value of previous step
                        IF (colIdx > rowIdx) THEN
                            m_sorX(i) = m_sorX(i) - (omega/m_matDiag(i)) * m_matOffDiag(j,i) * xx(lcc(j,i))
                        ELSEIF (colIdx < rowIdx) THEN
                            m_sorX(i) = m_sorX(i) - (omega/m_matDiag(i)) * m_matOffDiag(j,i) * m_sorX(lcc(j,i))
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO

        !...update xx
        DO i = 1,m_nOwnedDofs
            xx(i) = m_sorX(i)
        ENDDO

        !...release pointing to target
        CALL VecRestoreArrayF90(x, xx, p_err);
        
    END SUBROUTINE MatSOR_amat


    !...this subroutine performs v = A * u where v and u are ghosted array (i.e. including ghost cells)
    !...A is the coefficient matrix in the equation for solving intensity 
    SUBROUTINE matvec_ghosted(v, u)
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(inout) :: u
        REAL(real_p), DIMENSION(:), INTENT(out) :: v

        !...initialize v will be carried out in partial_matvec

        !...send values of u at owned cells to others who need as ghost, receive values of u at ghost cells from their owners
        CALL ghost_receive_begin(u)

        !...while waiting for communication, perform matvec on independent cells
        CALL partial_matvec(v, u, m_indCells, m_nIndCells)

        !...finishing communication
        CALL ghost_receive_end(u)
        
        !...once communication is completed, perform matvec on dependent elements
        CALL partial_matvec(v, u, m_depCells, m_nDepCells)

    END SUBROUTINE matvec_ghosted

    !...subroutine to perform v = K*u on a partial number of cells, used for ind/dep cells separately
    SUBROUTINE partial_matvec(v, u, cellList, listSize)
        IMPLICIT NONE
        REAL(real_p), DIMENSION(:), INTENT(in) :: u
        REAL(real_p), DIMENSION(:), INTENT(out) :: v
        INTEGER(int_p), DIMENSION(:), INTENT(in) :: cellList    !...list of cells to be performed v = K * u
        INTEGER(int_p), INTENT(in) :: listSize                  !...size of cellList

        !...local variables
        INTEGER(int_p) :: i, j, k

        !...loop over cell list
        DO k = 1,listSize
            i = cellList(k) !...cell (local) id
            v(i) = m_matDiag(i) * u(i) !...diagonal term, this also initializes
            !...loop over off-diagonal terms
            DO j = 1,nfcell(i)
                !...only multiply with off-diagonal if face j is an iterior face
                IF (bface(lcf(i,j)) .eq. 0) THEN
                    v(i) = v(i) + m_matOffDiag(j,i) * u(lcc(j,i))
                ENDIF
            ENDDO
        ENDDO
    END SUBROUTINE partial_matvec


    !...subroutine to compute each (local) row of the global matrix and the rhs vector
    SUBROUTINE compute_aMat_matrix_n_rhs
        IMPLICIT NONE

        !...local variables
        INTEGER(int_p) :: i, j, mat, f, currf, ibface, bcon, pol, sout, mat1
        INTEGER(int_p) :: ic_1, ic_2, pol1, pol2

        REAL(real_p) :: t_12, r_12, t_21, r_21, jface1, jface2, jface11, jface22
        REAL(real_p) :: velener_de, velener_r, n_omega1, n_omega2, disfun1

        REAL(real_p) :: temper, vol, gamma, sdotn, insdotn, vomega, vel, alpha, &
            rmuval, rxival, retval, inrmuval, inrxival, inretval, tempwall, &
            actin, actout, beta, mfp, jout, temper1

        pol = polar_type(iband_ctx, p_ctx)

        vel = gpvel(iband_ctx, pol)
        vomega = omega(si_ctx)
        rmuval = rmu(si_ctx)
        rxival = rxi(si_ctx)
        retval = ret(si_ctx)
        inrmuval = inrmu(si_ctx)
        inrxival = inrxi(si_ctx)
        inretval = inret(si_ctx)

        !...initialize, only need for diagonal terms
        m_matDiag(:) = zero
        m_matOffDiag(:,:) = zero
        m_rhs(:) = zero
        !...loop over cells in the list, for each cell i compute v(i) = K_ij * u_j where j = 1,nfaces(i)
        DO i = 1, m_nOwnedDofs
            mat = cell_type(i)          !! material_type
            !IF (mat == ALUMINUM) CYCLE
            temper = tnot(i)

            !f = checkpolar(mat, iband_ctx, p_ctx)
            !IF (f == 0) GO TO 100

            vol = volcell(i) * vomega

            CALL relaxtime(pol, mat, iband_ctx, temper, beta) !... what beta ?
            ! beta = 1.0D-6
            mfp = vel/beta

            !...global index for petsc vector m_rhs (even for CONJ case: ic_1 = i or ic_2 = i)
            m_idxs(i) = lc_to_pgc(i) - 1
            
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
                ! sdotn = sdotn_array(currf, si_ctx)
                ! insdotn = insdotn_array(currf, si_ctx)
                
                IF (bface(currf) == 0) THEN 
                    !...Interior face, no contribution to m_rhs
                    alpha = normdir(i,j)
                    sdotn = alpha * sdotn
                    insdotn = alpha*insdotn

                    m_matDiag(i) = m_matDiag(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))
                    m_matOffDiag(j,i) = -Max(0.0,-sdotn)*gamma*((insdotn/(sdotn+tiny)))

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
                        CALL reflec(si_ctx,currf,sdotn,sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)

                        !...store matrix
                        m_matDiag(i) = m_matDiag(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !...store m_rhs vector
                        m_rhs(i) = m_rhs(i) + dos*intensity(i,sout,icount_ctx,pol) * gamma * &
                                Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) &
                                +((one-dos)*jout) * gamma * Max(0.0,-sdotn)* &
                            ((insdotn/(sdotn+tiny)))/pi

                    CASE(ISOT)
                        !...store matrix
                        m_matDiag(i) = m_matDiag(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... store m_rhs vector
                        m_rhs(i) = m_rhs(i) + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    CASE (HFLUX) !!! assuming a virtual block of aluminium!!
                        !...store matrix
                        m_matDiag(i) = m_matDiag(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... store m_rhs vector
                        m_rhs(i) = m_rhs(i) + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))

                    CASE (ALSI)
                        !...store matrix
                        m_matDiag(i) = m_matDiag(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... store m_rhs vector
                        m_rhs(i) = m_rhs(i) + actin * gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))
                    
                    !!! symmetry boundary
                    CASE(SYMM)
                        !...store matrix
                        jout = jfacepos(ibface, pol, iband_ctx)
                        CALL reflec(si_ctx,currf,sdotn,sout)           !CALL reflec(si_ctx,ibface,sdotn,sout)
                        
                        m_matDiag(i) = m_matDiag(i) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                        !... store m_rhs vector
                        m_rhs(i) = m_rhs(i) + intensity(i,sout,icount_ctx,pol) * gamma * &
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
                    !         m_matDiag(ic_1) = m_matDiag(ic_1) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                    !         !... store m_rhs vector
                    !         m_rhs(ic_1) = m_rhs(ic_1) + velener_r*((jface1*r_12)/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) + &
                    !                     (jface22*t_21/pi)*gamma*Max(0.0,-sdotn)*((insdotn/(sdotn+tiny)))*velener_r

                    !     ELSE
                    !         !...store matrix
                    !         m_matDiag(ic_2) = m_matDiag(ic_2) + Max(0.0,sdotn)*gamma*((insdotn/(sdotn+tiny)))

                    !         !... store m_rhs vector
                    !         m_rhs(ic_2) = m_rhs(ic_2) + velener_r*((jface2*r_21)/pi)*gamma * Max(0.0,-sdotn)*((insdotn/(sdotn+tiny))) + &
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
            !f = checkpolar(mat1,iband_ctx, p_ctx)
            !IF (f /= 0) THEN
                CALL relaxtime(pol, mat1, iband_ctx, temper1, beta)
                !beta = 1.0D-6
                CALL iwalls(pol, mat1, iband_ctx, temper1, actout)

                !...store matrix
                m_matDiag(i) = m_matDiag(i) + (idt/beta + one)

                !...store m_rhs vector
                m_rhs(i) = m_rhs(i) + actout + idt *intensityone(i,si_ctx,icount_ctx,pol)/beta
            !ENDIF
        ENDDO
    END SUBROUTINE compute_aMat_matrix_n_rhs


    SUBROUTINE modify_aMat_rhs_for_dx_solving
        IMPLICIT NONE

        !...local variables
        INTEGER(int_p) :: i

        !...put intensity of previous iteration into m_u (we used m_u for this subroutine to avoid allocation of new memory)
        DO i = 1,m_nOwnedDofs
            m_u(i) = intensity(i, si_ctx, icount_ctx, p_ctx)
        ENDDO

        !...peform m_v = A * m_u where A is the matrix of the system of equations of intensity of current iteration
        CALL matvec_ghosted(m_v, m_u)

        !...subtract m_v from the m_rhs
        DO i = 1,m_nOwnedDofs
            m_rhs(i) = m_rhs(i) - m_v(i)
        ENDDO
    END SUBROUTINE modify_aMat_rhs_for_dx_solving


    SUBROUTINE compute_residual(sol, res)
        IMPLICIT NONE
        Vec, INTENT(IN) :: sol              !...solution (petsc) vector
        REAL(real_p), INTENT(OUT) :: res    !...L2 norm of residual vector

        !...local variables
        PetscScalar, Pointer, DIMENSION(:) :: ssol
        PetscErrorCode p_err
        INTEGER(int_p) :: i
        REAL(real_p) :: temp1, temp2, temp1_sum, temp2_sum
        INTEGER(int_p), PARAMETER :: root_proc = 0

        CALL VecGetArrayReadF90(sol, ssol, p_err)

        DO i = 1,m_nOwnedDofs
            m_u(i) = ssol(i)
        ENDDO

        !... m_v = A * m_u
        CALL matvec_ghosted(m_v, m_u)

        !...subtract rhs from m_v
        DO i = 1,m_nOwnedDofs
            m_v(i) = m_v(i) - m_rhs(i)
        ENDDO

        !...compute L2 norm of m_v and m_rhs
        temp1 = 0.d0
        temp2 = 0.d0
        DO i = 1,m_nOwnedDofs
            temp1 = temp1 + m_v(i) * m_v(i)
            temp2 = temp2 + m_rhs(i) * m_rhs(i)
        ENDDO

        !...reduction sums to root
        CALL MPI_REDUCE(temp1, temp1_sum, 1, MPI_REAL8, MPI_SUM, root_proc, m_comm, p_err)
        CALL MPI_REDUCE(temp2, temp2_sum, 1, MPI_REAL8, MPI_SUM, root_proc, m_comm, p_err)

        !... root computes norm of residual
        IF (m_rank .eq. root_proc) THEN
            temp1_sum = SQRT(temp1_sum)
            temp2_sum = SQRT(temp2_sum)
            res = temp1_sum/temp2_sum
        ENDIf

        !... broadcast res computed by root to other ranks
        CALL MPI_BCAST(res, 1, MPI_REAL8, root_proc, m_comm, p_err)

    END SUBROUTINE compute_residual


    !...set aMat rhs to PETSc vector
    SUBROUTINE set_petsc_rhs
        IMPLICIT NONE

        CALL VecSetValues(mp_rhs, m_nOwnedDofs, m_idxs, m_rhs, INSERT_VALUES, mp_ierr)

        CALL VecAssemblyBegin(mp_rhs, mp_ierr)
        CALL VecAssemblyEnd(mp_rhs, mp_ierr)
    END SUBROUTINE set_petsc_rhs


    !...set s = intensity of previous iteration, for solve_option = 1/3 of aMat
    SUBROUTINE set_sol_prev_itsy(s)
        IMPLICIT NONE
        Vec s

        !...internal variables
        INTEGER(int_p) :: i
        PetscScalar, Pointer, DIMENSION(:) :: ss
        PetscErrorCode p_err

        !...point to petsc vector allowing ss to read/write to s
        CALL VecGetArrayF90(s, ss, p_err); CHKERRA(p_err)
        DO i = 1,m_nOwnedDofs
            ss(i) = intensity(i, si_ctx, icount_ctx, p_ctx)
        ENDDO

        !...release pointing target
        CALL VecRestoreArrayF90(s, ss, p_err); CHKERRA(p_err)
        
    END SUBROUTINE set_sol_prev_itsy

END MODULE amat_noBatched
