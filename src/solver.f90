PROGRAM Main

   USE PRECISIONS
   USE VARIABLES
   USE GRID
   USE CONSTANTS, only : one, four, pi, zero, two, notsotiny, tiny, k_si

   !...modules corresponding to method used for solving intensity
   USE petsc_matAssembled
   USE amat_noBatched
   USE amat_batched

   IMPLICIT NONE

   INTEGER(int_p) :: i, iband, pol, tindex, iter, j , mat, f, ic, ifc, ic_1, ic_2, ibface, p
   INTEGER(int_p) :: cell
   REAL(real_p) :: phiguess, gval, flux_net

   LOGICAL :: exists = .true.

   INTEGER(int_p) :: icount
   INTEGER :: id, flag
   INTEGER, PARAMETER :: root_proc = 0

   REAL(real_p) :: t_time

   CHARACTER(256) :: model_setup_bte_filename, user_polar_filename
   CHARACTER(8) :: arg3, arg4, arg5

   INTEGER, DIMENSION(0:1) :: dims
   LOGICAL, DIMENSION(0:1) :: periods = (/.FALSE., .FALSE./)
   LOGICAL, DIMENSION(0:1) :: remain_dims
   LOGICAL :: reorder = .TRUE.
   INTEGER(int_p) :: nProcsRow, nProcsCol
   INTEGER :: request

   REAL(real_p) :: avg_ttop
   REAL(real_p) :: timeStep_max, timeStep_min
   !REAL(real_p) :: disord_max, disord_min
   REAL(real_p) :: previousTimer
   REAL(real_p) :: m_mvTimer_max, m_mvTimer_min, m_precondTimer_max, m_precondTimer_min

   !...runtime parameters: input files, n procs for band-based, n procs for cell-based, matrix method, batch size
   CALL getarg(1, model_setup_bte_filename)
   CALL getarg(2, user_polar_filename)
   CALL getarg(3, arg3) !...np for band-based
   CALL getarg(4, arg4) !...np for cell-based
#ifdef USE_AMAT_BATCHED
   CALL getarg(5, arg5) !...batch size
#endif
   !CALL getarg(6, arg6)

   !...convert char to numeric values
   READ(arg3, *) nProcsRow !...n processes for band-based
   READ(arg4, *) nProcsCol !...n processes for cell-based
#ifdef USE_AMAT_BATCHED
   READ(arg5, *) batchSize
#endif

   !...initialize PETSc
   CALL PetscInitialize(PETSC_NULL_CHARACTER, mp_ierr)
   IF (mp_ierr .ne. 0) THEN
      WRITE(*,'(a)') 'Unable to initialize PETSc!'
      STOP
   ENDIF

   CALL MPI_COMM_SIZE(PETSC_COMM_WORLD, numprocs, mp_ierr)
   CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, mp_ierr)
   IF ((nProcsRow * nProcsCol) .ne. numprocs) THEN
      WRITE(*,'(a)') 'Number of rows * number of cols must be equal to np!'
      CALL MPI_Abort(PETSC_COMM_WORLD, -1, mpierr)
   ENDIF

   !...create Cartesian communicator
   dims(0) = nProcsRow
   dims(1) = nProcsCol
   CALL MPI_Cart_create(PETSC_COMM_WORLD, 2, dims, periods, reorder, comm_cart, mp_ierr)

   !...size and rank of Cartesian is the same as mpi_comm_world
   CALL MPI_COMM_RANK(comm_cart, rank_cart, mp_ierr)
   CALL MPI_COMM_SIZE(comm_cart, numprocs_cart, mp_ierr)

   !...sub-grid communicator for cell-based (processes in one row have the same comm_cell)
   remain_dims(0) = .FALSE.
   remain_dims(1) = .TRUE.
   CALL MPI_Cart_sub(comm_cart, remain_dims, comm_cell, mp_ierr)
   CALL MPI_Comm_size(comm_cell, numprocs_cell, mp_ierr)
   CALL MPI_Comm_rank(comm_cell, rank_cell, mp_ierr)
   IF (numprocs_cell .ne. nProcsCol) THEN
      WRITE(*,'(a)') 'Error: n of processes in column direction is not equal to n processes used for cell-based parallelization!'
      STOP
   ENDIF

   !...sub-grid communicator for band-based (processes in one column have the same comm_band)
   remain_dims(0) = .TRUE.
   remain_dims(1) = .FALSE.
   CALL MPI_Cart_sub(comm_cart, remain_dims, comm_band, mp_ierr)
   CALL MPI_Comm_size(comm_band, numprocs_band, mp_ierr)
   CALL MPI_Comm_rank(comm_band, rank_band, mp_ierr)
   IF (numprocs_band .ne. nProcsRow) THEN
      WRITE(*,'(a)') 'Error: n of processes in row direction is not equal to n processes used for band-based parallelization!'
      STOP
   ENDIF

   !...file for writing total time
   IF (rank_cart .eq. root_proc) THEN
      IF (debug_level > 0) THEN
         OPEN(UNIT=io_debug, FILE="debug.out", status="unknown")
         WRITE(io_debug, '(a)') "Starting solver (main)"
      ENDIF

      !...display model parameters
      WRITE(*, '(a)') '=================================================='
      WRITE(*, '(a,a)') 'Model file: ', model_setup_bte_filename
      WRITE(*, '(a,a)') 'Parameter file: ', user_polar_filename
      WRITE(*, '(a,i0,a,i0)') 'N processes for band-based = ',nProcsRow,', N processes for cell-based = ',nProcsCol
#ifdef USE_PETSC_MATASSEMBLED
      WRITE(*,'(a)') 'PETSc matrix-assembled method is used for intensity solving'
#elif USE_AMAT_NOBATCHED
      WRITE(*,'(a)') '(no batching) aMat method is used for intensity solving'
#elif USE_AMAT_BATCHED
      WRITE(*,'(a)') 'aMat method, with batching, is used for intensity solving'
      WRITE(*,'(a,i0)') 'batch size = ',batchSize
#else
      WRITE(*,'(a)') 'Incorrect solving method is chosen for compiling!'
      STOP
#endif
      
      WRITE(*, '(a)') '=================================================='
   ENDIF

   WRITE(nProcsName, '(i0)') numprocs  !...total processes used including processes for band-based in the combined method
   WRITE(rankName, '(i0)') rank        !...the postprocess is called only by root of rank_band, so only processes having rank-band=0 will do postprocess

   CALL profiler_init(timeStep_timer)

   !...reads model parameters in user_polar_filename
   IF (rank_cart .eq. root_proc) THEN
      INQUIRE(file=user_polar_filename, EXIST = exists)
      IF (.not. exists) THEN
         WRITE(*,'(a)') "Parameter file could not be found!"
         STOP
      ENDIF
      
      OPEN(UNIT=io_user_input, FILE=user_polar_filename, status="old")
      READ(io_user_input,*) initial_temp
      READ(io_user_input,*) laserpower    !...not be broadcasted
      READ(io_user_input,*) layer         !...not be broadcasted
      READ(io_user_input,*) spotsize      !...not be broadcasted
      READ(io_user_input,*) probe_size    !...for ALSI boundary condition
      READ(io_user_input,*) coldwall
      READ(io_user_input,*) hotwall
      READ(io_user_input,*) nbands        !! change input file
      READ(io_user_input,*) dos
      READ(io_user_input,*) nphi
      READ(io_user_input,*) ntheta
      READ(io_user_input,*) dt_1
      READ(io_user_input,*) tmax
      READ(io_user_input,*) timestep_period
      READ(io_user_input,*) n
      READ(io_user_input,*) time_pulse
      READ(io_user_input,*) mod_freq
      READ(io_user_input,*) tol_outer
      READ(io_user_input,*) rel           !...not be broadcasted

      CLOSE(UNIT=io_user_input)
      WRITE(*,'(a,e15.5)') 'Outer tolerance = ', tol_outer
      WRITE(*,'(a,i0)') "Total processes= ", numprocs
      
      ndir = Int(ntheta * nphi) !from dircalc
      dphi = two * pi/nphi  ! dircalc
      dtheta = pi/ntheta  !dircalc
   ENDIF

   !... root of Cartesian broadcasts model parameters to all other ranks
   CALL MPI_BCAST(initial_temp, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(coldwall, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(hotwall, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(nbands,1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(dos, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(nphi, 1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(ntheta, 1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(dt_1, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(tmax, 1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(timestep_period, 1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(n, 1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(time_pulse, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(mod_freq, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(tol_outer, 1, MPI_REAL8, root_proc, comm_cart, mpierr)

   !...ndir, dphi, dtheta are computed by p0
   CALL MPI_BCAST(ndir, 1, MPI_INTEGER, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(dphi, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(dtheta, 1, MPI_REAL8, root_proc, comm_cart, mpierr)

   CALL MPI_BCAST(laserpower, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(spotsize, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(rel, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(layer, 1, MPI_REAL8, root_proc, comm_cart, mpierr)
   CALL MPI_BCAST(probe_size, 1, MPI_REAL8, root_proc, comm_cart, mpierr)

   !...all ranks read global mesh
   CALL grid_and_boundary(model_setup_bte_filename)

   IF (rank_cart .eq. root_proc) THEN
      WRITE(*,'(a)') "Finish grid_and_boundary"
      IF (debug_level > 0) WRITE(io_debug,'(a)') "Finish grid_and_boundary"
   ENDIF

   !...for parallel in bands
   ALLOCATE(displacement(0:numprocs_band - 1))
   ALLOCATE(recievecount(0:numprocs_band - 1))
   ALLOCATE(mpistatus(MPI_STATUS_SIZE))
   ALLOCATE(mpirequest(MPI_REQUEST_NULL))

   lowband = 1 + (rank_band * nbands)/numprocs_band     ! indexing starts from 1 in fortran
   highband = ((rank_band + 1) * nbands)/numprocs_band  ! num process should be multiple of bands!

   !...allocate and initialize variables, also create PETSc matrix, rhs and sol vectors
   CALL Initialize

   IF (rank_cart .eq. root_proc) THEN
      WRITE(*,'(a)')'Finish initialize'
      IF (debug_level > 0) WRITE(io_debug,'(a)') "Finish initialize"
   ENDIF

   !...initialize control angle discrete coordinate (in snweights.f90)
   CALL dircalc

   !...set Gauss points and weights xi() and wi() for 20-point formula
   CALL gauss_points

   CALL polarization

   !...compute band averaged velocity of each frequency band
   CALL avgvel

   !...overwrite temperature on boundary face (read in from input) by initial_temp
   DO ibface = 1, nbcfaces
      IF(bctype(ibface) .ne. ISOT) temp_bc(ibface) = initial_temp
   ENDDO

   temp_bcone(:) = temp_bc(:)

   DO i = 1, ncells
      tnot(i) = initial_temp
   ENDDO

   DO iband = 1,nbands
      !... np_max is max of polarization, computed in grid_and_boundary
      DO j = 1, np_max
         pol = polar_type(iband,j)
         IF (pol .eq. 0) CYCLE
         DO i = 1,ncells
            mat = cell_type(i)
            !IF(mat == ALUMINUM) CYCLE
            !f = checkpolar(mat,iband, j)
            !IF (f == 0) CYCLE
            !... compute intensity from temperature
            CALL iwalls(pol, mat, iband, tnot(i), phiguess)
            
            gval =  four*pi*phiguess   !...4*pi is the total direction angle
            gc(i, j, iband) = gval     !...integrated intensity over all directions: 3.18
            gcone(i, j, iband) = gval
         ENDDO
      ENDDO
   ENDDO

   DO iband = lowband, highband
      icount = iband - lowband + 1
      DO j = 1, np_max
         pol = polar_type(iband,j) !! pol = {0, 1, 2, 3, 4}
         IF (pol .eq. 0) CYCLE
         DO i = 1,ncells
            mat = cell_type(i)
            !f = checkpolar (mat,iband, j)
            !IF (f == 0) CYCLE
            CALL iwalls(pol, mat, iband, tnot(i), phiguess)
            intensity(i, :, icount, j) = phiguess
            intensityone(i, :, icount, j) = phiguess
         ENDDO
      ENDDO
   ENDDO

   !CALL variabledt !...comment out on Jan 30 according to Siddharth's implicit_memory code
   sep = 1 !...according to Siddharth's implicit_memory code

   npulse = zero

   !CALL profiler_start(timeStep_timer)
   DO tindex = 1,tmax  ! Time Stepping
      previousTimer = timeStep_timer%seconds
      CALL profiler_start(timeStep_timer)
      ntimep = MOD((tindex-1), timestep_period)
      dt = dt_1 * sep**(ntimep)
      idt = one/dt
      IF (ntimep == 0) THEN
         npulse = npulse + 1
      ENDIF
      IF (sep == 1) THEN
         t_time = ((npulse-1) * time_pulse) + (ntimep*dt_1)
         IF (rank_cart .eq. root_proc) WRITE(*,'(f20.15)') t_time
      ELSE
         t_time = ((npulse-1) * time_pulse) + ((1 - (sep**(ntimep+1)))/(1-sep))*dt_1
      ENDIF

      !! find the power (amount of heat putting into a laser)
      IF (rank_band .eq. root_proc) THEN 
         residual = one
         qdblprime = laserpower * (1 + sin(mod_freq * t_time * two * pi))
         !!!!!! silicon thermal conductivity
         !IF (rank_cell .eq. root_proc) WRITE (*,'(a,i0,1x,f20.15)') "tindex= ", tindex, qdblprime
         !WRITE (*,*) "tindex= ", rank_band, rank_cell, tindex, qdblprime, t_time, mod_freq, laserpower
      ENDIF

      !... k_sil is conductivity of silicon
      k_sil(:) = k_si

      !...flag to indicate if convergence (of outer loop) occurs, all ranks set it to NO
      flag = 0

      !...sum of residual which is reduced to rank 0
      residual_sum = 0.0

      !...start outer loop
      DO iter = 1,max_iter
         !... compute jfacepos (another type of integrated intensity)
         CALL intintesity
         
#ifdef USE_PETSC_MATASSEMBLED
         CALL disord ! CADOM for intensity of all bands and directions!
#elif USE_AMAT_NOBATCHED
         CALL disord_aMat(iter)
#elif USE_AMAT_BATCHED
         CALL disord_aMat_bch
#endif

         !!! NEED TO GATHER GC for residual and detT calcualation
         displacement(0) = 0
         recievecount(0) = (highband - lowband + 1) * ncells * np_max
         DO id = 1, (numprocs_band - 1)
            recievecount(id)= (((id + 1) * nbands)/numprocs_band - &
               (id * nbands)/numprocs_band) * ncells * np_max ! if no of bands divisible by numpros, then nband/numporces
            displacement(id) = displacement(id-1) + recievecount(id-1)
         ENDDO

         IF (rank_band .eq. root_proc) THEN
            CALL MPI_GATHERV(MPI_IN_PLACE, ((highband - lowband + 1) * ncells * np_max), MPI_REAL8, &
               gc(1,1,lowband), recievecount, displacement, MPI_REAL8, root_proc, comm_band, mp_ierr)
         ELSE
            CALL MPI_GATHERV(gc(1,1,lowband), ((highband - lowband + 1) * ncells * np_max), MPI_REAL8, &
               gc(1,1,lowband), recievecount, displacement, MPI_REAL8, root_proc, comm_band, mp_ierr)
         ENDIF

         CALL intintesity !... integrated intensity

         !!gather flux for use of postprocess
         displacement(0)=0
         recievecount(0)=(highband - lowband + 1) * nbcfaces * np_max
         DO id = 1, (numprocs_band - 1)
            recievecount(id)= (((id + 1) * nbands)/numprocs_band - &
               (id * nbands)/numprocs_band) * nbcfaces * np_max ! if no of bands divisible by numpros, then nband/numporces
            displacement(id) = displacement(id-1) + recievecount(id-1)
         ENDDO

         CALL MPI_GATHERV(jfacepos_s(1,1,lowband), ((highband - lowband + 1) * nbcfaces * np_max), MPI_REAL8, &
            jfacepos(1,1,lowband), recievecount, displacement, MPI_REAL8, root_proc, comm_band, mpierr)

         CALL MPI_GATHERV(jfaceneg_s(1,1,lowband), ((highband - lowband + 1) * nbcfaces * np_max), MPI_REAL8, &
            jfaceneg(1,1,lowband), recievecount, displacement, MPI_REAL8, root_proc, comm_band, mpierr)

         ! MPI GATHER TO GET flux(:,:,:) to root_process
         CALL MPI_GATHERV(flux_s(1,1,lowband), ((highband - lowband + 1) * nbcfaces * np_max), MPI_REAL8, &
            flux(1,1,lowband), recievecount, displacement, MPI_REAL8, root_proc, comm_band, mpierr)
   
         !...after gather to root_proc of comm_band (there are numprocs_cell processes that are root in comm_band), 
         !...now processes that are root of comm_band compute the temperature for the cells they own
         IF (rank_band .eq. root_proc) THEN
            CALL tempcalc ! calculate temperature (First Law)

            !!!!!!!!!!!!!!writing iteration with tempcalc
            IF (tmax .lt. 15) THEN
               DO i = 1,nbcfaces
                  ifc = bf_to_f(i)
                  ic = lfc (ifc,1)
                  flux_net = zero
                  Do iband = 1, nbands
                     DO p = 1, np_max
                        flux_net = flux_net - flux(i,p,iband)
                     ENDDO
                  ENDDO
                  ! IF(bctype(i) == ALSI) THEN
                  !    WRITE(61,204) tindex, iter,xf(ifc),yf(ifc),ttop(i),temp_bc(i), tnot(ic), flux_net, qdblprime
                  ! ENDIF
               ENDDO
            ENDIF

            !...reduce sum of residual (before taking square root) across ranks of comm_cell to the root (of comm_cell)
            CALL MPI_REDUCE(residual, residual_sum, 1, MPI_REAL8, MPI_SUM, root_proc, comm_cell, mp_ierr)

            !...only root (of comm_cell, thus this one is the root of both comm_band and comm_cell) takes square root of sum of residual
            IF (rank_cell .eq. root_proc) THEN
               residual_sum = SQRT(residual_sum)
               WRITE(*,'(a,i0,a,i0,a,e15.8)') "Time step ",tindex,", iteration ",iter,", residual= ",residual_sum
               !WRITE(*,'(a,i0,a,i0,a,e15.8,a,i0)') "Time step ",tindex,", iteration ",iter,", residual= ",residual_sum, ", total iterations of linear solving= ", m_nIts_total
               !...since the root get the sum reduction of residual, it checks the convergence condition
               IF (iter > min_iter) THEN
                  IF(residual_sum < tol_outer)THEN
                     flag = 1
                     !WRITE(*,'(a,i0,a,e15.8)') "Finished time step ",tindex,", residual= ",residual_sum
                  ENDIF
               ENDIF
               IF (iter == max_iter) THEN
                  WRITE(*,'(a,i0)') "Number of (outer) iterations exceeded max_iter (i.e. did not converge) for time step ",tindex
               ENDIF
            ENDIF

            !...broadcast the flag from the root (of comm_cell) to all processes to stop iteration if flag = 1
            !...we may want to combine this with the broadcast in comm_band into one broadcast in comm_cart
            CALL MPI_BCAST(flag, 1, MPI_INTEGER, root_proc, comm_cell, mpierr)
         ENDIF

         !...now broadcast the flag from all roots (of comm_band) to all processes in comm_band
         CALL MPI_BCAST(flag, 1, MPI_INTEGER, root_proc, comm_band, mpierr)

         ! tnot broadcast for update of variables!!
         CALL MPI_BCAST (tnot(1), ncells, MPI_REAL8, root_proc, comm_band, mpierr)

         ! temp_bc broadcast for update of variables!!
         CALL MPI_BCAST (temp_bc(1), nbcfaces, MPI_REAL8, root_proc, comm_band, mpierr)

         !...Nov 10, 2022: Siddharth moves this exit after the broadcast of tnot() and temp_bc() to fix the discrepancies between sequential and band-based parallel after time step= 1
         IF (flag .eq. 1) THEN
            EXIT
         ENDIF
      ENDDO

      ! Updates Solution for next time-step
      CALL update

      displacement(0) = 0
      recievecount(0) = (highband - lowband + 1) * np_max
      DO id = 1, (numprocs_band - 1)
         recievecount(id) = (((id+1)*nbands)/numprocs_band - (id*nbands)/numprocs_band)*np_max ! if no of bands divisible by numpros, then nband/numporces
         displacement(id) = displacement(id-1) + recievecount(id-1)
      ENDDO

      CALL profiler_stop(timeStep_timer) !...exclude the timeStep_timer from post-processing time
      IF (rank_cart .eq. root_proc) THEN
         WRITE(*,'(a,i0,a,f15.5,a,i0,a,f15.5,a)') '(root process) time step ', tindex, ' took ', (timeStep_timer%seconds - previousTimer)/3600, '(h); n iterations= ', iter, &
            ', accumulated time = ', timeStep_timer%seconds/3600, ' (h)'
      ENDIF 

      IF (rank_band .eq. root_proc) THEN
         IF ((mod(tindex,freq_time_output)==0).or.(tindex == 1)) then  
            CALL postprocess(tindex)
         ENDIF

         !...each process (who has rank_bank = 0) computes its local sums of area-weighted ttop and area
         !...then mpi_reduce sum to process who has rank_cell = 0, then this process compute the average ttop
         !CALL avg_ttop_calc(avg_ttop)
         !IF (rank_cell .eq. 0) WRITE(*,*) 'Averaged ttop = ', avg_ttop

         DO i = 1,nbcfaces
            IF((bctype(i) == ALSI) .or. (bctype(i)== HFLUX))THEN
               ifc = bf_to_f(i)
               ic = lfc (ifc,1)
               ! IF (cell_type(ic) == ALUMINUM) ic = lfc (ifc,2)
               !qnot = qdblprime*exp(-(2*yf(ifc)*yf(ifc))/(spotsize*spotsize)) !...TODO: WE NEED TO CHANGE THIS (ONLY APPLICABLE TO 2D)
               !IF (yf(ifc) .gt. spotsize) qnot = zero 
               flux_net = zero    
               Do iband = 1, nbands
                  DO p = 1, np_max
                     flux_net = flux_net + flux(i,p,iband)
                  ENDDO
               ENDDO
               !...this is according to Siddharth's implicit_memory code (need to confirm w/ Siddharth b/c sep=1 is hard code, which makes different with before)
               IF (sep == 1) THEN
                  t_time=((npulse-1) *time_pulse)+(ntimep*dt_1)
               ELSE
                  t_time = ((npulse-1) *time_pulse) + ((1 - (sep**(ntimep+1)))/(1-sep))*dt_1 !...this is used to be
               ENDIF
            ENDIF      
         ENDDO
      ENDIF
   ENDDO ! End time loop

   !...reduction to get max/min timers among processes
   CALL MPI_REDUCE(timeStep_timer%seconds, timeStep_max, 1, MPI_REAL8, MPI_MAX, root_proc, comm_cart, mpierr)
   CALL MPI_REDUCE(timeStep_timer%seconds, timeStep_min, 1, MPI_REAL8, MPI_MIN, root_proc, comm_cart, mpierr)

   IF (rank_cart .eq. root_proc) THEN
      !...timers measured by main():
      WRITE(*,'(a,f15.5,a,f15.5,a,f15.5)') 'max/min of complete total time (not include post-process) = ', timeStep_max/3600, ' (h) /', &
         timeStep_min/3600,' (h), ratio = ', timeStep_max/timeStep_min
   ENDIF

   IF (rank_cart .eq. root_proc) THEN
      IF (debug_level > 0) THEN
         WRITE(io_debug,'(a)') "Finishing solver (main)" 
         CLOSE(UNIT=io_debug)
      ENDIF
   ENDIF

   201 FORMAT (I6,E16.8,E12.5,E12.5,F16.8,F16.8, F16.8,F16.8,F16.8,E16.8)
   202 FORMAT (I6,E15.8,E12.3,F16.8,E15.7,E15.7)
   203 FORMAT (I6,E15.8,E12.5,E12.5,E12.5,F16.8,F16.8,F16.8, F16.8,F16.8,E16.8)
   204 FORMAT (I6,I6,E15.5,E15.5, F18.12,F18.12,F18.12,E16.8, E16.8)
   205 FORMAT (I6, E16.8,E16.8,E16.8, E16.8,E16.8,F16.8,F16.8 )

   !... free PETSc vectors
   CALL destroy_petsc_rhs_n_sol

   !... free PETSc matrix
#ifdef USE_PETSC_MATASSEMBLED
   CALL destroy_petsc_mat
#elif USE_AMAT_NOBATCHED
   CALL finalize_aMat
#elif USE_AMAT_BATCHED
   CALL finalize_aMat_bch
#endif

   !...free PETSc KSP
   CALL destroy_petsc_ksp

   !CALL MPI_FINALIZE (mpierr) !... this is done inside PetscFinalize()
   CALL PetscFinalize(mp_ierr)

END PROGRAM Main