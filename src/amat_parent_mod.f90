! \brief module to solve using assembly-free aMat (parent of amat_noBatched_mod and amat_batched_mod)
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! Solve the linear systems of intensity using assembly-free method.
! Parent module of non-batching and batching modules
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
MODULE amat_parent
    USE petsc_solver
    USE profiler
    IMPLICIT NONE

    !...scatter map, global variables for anyone who uses amat_parent
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendCellIds      !... local Ids of cells that I need to send to other ranks
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendCellCounts   !... number of cells that I need to send to each rank
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendCellOffset   !... offset of sendCellCounts
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_recvCellCounts   !... number of cells that I will receive from each rank
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_recvCellOffset   !... offset of recvCellCounts
    INTEGER(int_p) :: m_totalSend                                   !... total number of cells that I need to send to other ranks
    INTEGER(int_p) :: m_totalRecv                                   !... total number of cells that I receive from other ranks (= number of ghost cells)

    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_sendBuf, m_recvBuf             !... buffer for sending and receiving
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_sendRankIds, m_recvRankIds   !... list of rank ids that I need to send to/ receive from
    INTEGER(int_p) :: m_nProcsSend, m_nProcsRecv                                !... number of ranks that I will send to / receive from
    
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_request          !... request of isend/irecv for use in mpi_wait; 2022_09_07: this is needed to be global because it is used in both ghost_receive_begin() and ghost_receive_end()
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: m_statusArray
    
    INTEGER(int_p) :: m_globalDofBegin, m_globalDofEnd              !... start and end of global cell Id owned by my rank, used to check if cell owed by my in building mp_matPC

    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_matDiag            !...diagonal terms of local row i of the global matrix, used in afm approach
    REAL(real_p), DIMENSION(:,:), ALLOCATABLE :: m_matOffDiag       !...off-diagonal terms of local row i
    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_rhs                !...rhs vector, added Oct 18 2021
    PetscInt, DIMENSION(:), ALLOCATABLE :: m_idxs                   !...global index for setting m_rhs to Petsc vector, added Oct 18 2021

    REAL(real_p), DIMENSION(:), ALLOCATABLE :: m_u, m_v             !...ghosted-include vectors v = K*u where K is the coefficient matrix

    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_IndCells         !...list of independent cells
    INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: m_depCells         !...list of dependent cells
    INTEGER(int_p) :: m_nIndCells, m_nDepCells                      !...number of independent/dependent cells
    !TYPE(profiler_t) :: m_mvTimer   !...timers of matvec
    !TYPE(profiler_t) :: m_precondTimer

    CONTAINS

    !...subroutine to identify independent and dependent cells
    SUBROUTINE identifyIndCells
        IMPLICIT NONE

        INTEGER(int_p) :: i, j, ierr
        INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: depCellFlag
        INTEGER(int_p) :: depCount, indCount

        ALLOCATE(depCellFlag(m_nOwnedDofs), stat=ierr)
        depCellFlag = 0

        !...search for dependent cells
        m_nDepCells = 0
        DO i = 1,m_nOwnedDofs
            DO j = 1,nfcell(i)
                !...if a neighbor of cell i is a ghost cell, then i is dependent
                IF (lcc(j,i) > m_nOwnedDofs) THEN
                    m_nDepCells = m_nDepCells + 1
                    depCellFlag(i) = 1
                    EXIT
                ENDIF
            ENDDO
        ENDDO

        !...number of independent cells
        m_nIndCells = m_nOwnedDofs - m_nDepCells

        !...allocate lists of dependent/independent cells
        IF (m_nIndCells > 0) ALLOCATE(m_indCells(m_nIndCells), stat=ierr)
        IF (m_nDepCells > 0) ALLOCATE(m_depCells(m_nDepCells), stat=ierr)

        !...specify dependent/independent cells
        depCount = 0
        indCount = 0
        DO i = 1,m_nOwnedDofs
            IF (depCellFlag(i) == 1) THEN
                depCount = depCount + 1
                m_depCells(depCount) = i
            ELSE
                indCount = indCount + 1
                m_indCells(indCount) = i
            ENDIF
        ENDDO

        !...check if the above calculation is correct
        IF (depCount .ne. m_nDepCells) THEN
            write(*,'(a,i0)') "identifyCells: error in computing dependent cells, nDepCells= ",m_nDepCells
            stop
        ENDIF

        !...free space
        IF (ALLOCATED(depCellFlag)) DEALLOCATE(depCellFlag)

    END SUBROUTINE identifyIndCells
    
END MODULE aMat_parent
