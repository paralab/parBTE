! \brief compute avarage top temperature
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! subroutine to compute average ttop for boundary faces within the probe area, applied for ALSI boundary condition
! this subroutine should be called only for processes having rank_band = root
SUBROUTINE avg_ttop_calc(avg_ttop)
   USE precisions
   USE variables, ONLY : weighted_ttop_local_sum, bfArea_local_sum, weighted_ttop_global_sum, &
      bfArea_global_sum, ttop, probe_size, rank_cell, comm_cell, mpierr
   USE grid, ONLY: areaf, bf_to_f, nbcfaces, xf, yf, zf, threed, bctype, ALSI
   USE MPI
   
   IMPLICIT NONE
   
   REAL(real_p), intent(out) :: avg_ttop
   INTEGER :: i, currf
   REAL(real_p) :: radius
   
   !...compute the local sum of weighted (by area) ttop, and boundary-face areas inside the probe
   weighted_ttop_local_sum = 0.0
   bfArea_local_sum = 0.0
   DO i = 1,nbcfaces
      IF (bctype(i) .eq. ALSI) THEN
         currf = bf_to_f(i)
         IF (threed) THEN
            radius = SQRT((xf(currf) * xf(currf)) + (yf(currf) * yf(currf)) + (zf(currf) * zf(currf)))
         ELSE
            radius = SQRT((xf(currf) * xf(currf)) + (yf(currf) * yf(currf)))
         ENDIF
         IF (radius .le. probe_size) THEN
            weighted_ttop_local_sum = weighted_ttop_local_sum + ttop(i) * areaf(currf)
            bfArea_local_sum = bfArea_local_sum + areaf(currf)
            !write(*,*)'(i,...)',i,ttop(i),areaf(currf),weighted_ttop_local_sum,bfArea_local_sum
         ENDIF
      ENDIF
   ENDDO

   !...mpi reduce to get global sums of weighted ttop and boundary-face area
   CALL MPI_REDUCE(weighted_ttop_local_sum, weighted_ttop_global_sum, 1, MPI_REAL8, MPI_SUM, 0, comm_cell, mpierr)
   CALL MPI_REDUCE(bfArea_local_sum, bfArea_global_sum, 1, MPI_REAL8, MPI_SUM, 0, comm_cell, mpierr)

   !...return average ttop
   IF (rank_cell .eq. 0) THEN
      avg_ttop = weighted_ttop_global_sum / bfArea_global_sum
      !write(*,*)'ttop=...',avg_ttop, weighted_ttop_global_sum, bfArea_global_sum
   ELSE
      avg_ttop = 0.0
   ENDIF

END SUBROUTINE avg_ttop_calc

