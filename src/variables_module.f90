#include <petsc/finclude/petscksp.h>
Module VARIABLES
  Use PRECISIONS
  USE profiler
  use petscksp
  Implicit None
  SAVE


  INTEGER(int_p) :: nphi, ntheta, tmax
  INTEGER(int_p) :: ndir, nbands
  INTEGER(int_p) :: nla, nta,nlo,nto ! No of bands in la and ta polarizations
  INTEGER (int_p) :: no_b_steps, np_max
  
  !...INTEGER(int_p), Allocatable ::  IA(:), JA(:)
  !...INTEGER(int_p), Allocatable, Dimension(:) :: jlu, ju, jr
  INTEGER(int_p), Allocatable, Dimension(:) :: polar
  
  INTEGER(int_p), ALLOCATABLE, DIMENSION (:) ::  no_polar, total_bands_in_step
  INTEGER(int_p), ALLOCATABLE, DIMENSION (:,:) :: polar_type !, mat_polar
  !INTEGER(int_p), ALLOCATABLE, DIMENSION (:,:,:) :: checkpolar
  
  INTEGER(int_p) :: batchSize !...number of systems Ax=b in one batch, get from command line
  INTEGER(int_p), ALLOCATABLE, DIMENSION(:) :: nBatchesPerBand !...number of batches (of systems Ax=b) per band
  INTEGER(int_p), ALLOCATABLE, DIMENSION(:,:) :: si_nonzeroSystemPerBand, p_nonZeroSystemPerBand

  !...for profiling
  TYPE(profiler_t) :: timeStep_timer !, disord_timer !...timers for complete time step, solving for intensity + flux + gc
  
  INTEGER (int_p) :: iter_cal, timeindex
  
  !!!!! for variable dt
  ! REAL(real_p) :: dt_1

  INTEGER (int_p) :: num_material

  !!! to get interior flux
  INTEGER (int_p) :: ref_cross_facefar, ref_trans_facefar,ref_cross_facenear, ref_trans_facenear
   
  !! MPI VARIABLES
  INTEGER :: lowband, highband
  INTEGER :: numprocs, rank, mpierr, stat
  INTEGER, ALLOCATABLE :: recievecount(:), displacement(:)  !array for MPI GATHER, to get all values in root process   
  INTEGER, ALLOCATABLE, DIMENSION (:) :: mpistatus, mpirequest
  INTEGER :: rank_cell, rank_band !...my rank in cell-based and band-based communicators, respectively
  INTEGER :: numprocs_cell, numprocs_band !...size of cell-based and band-based communicators, respectively
  INTEGER :: comm_cell, comm_band !...communicators of cell-based and band-based parallelization, respectively
  INTEGER :: comm_cart, rank_cart, numprocs_cart

  ! File numbers used throughout the code
  INTEGER(int_p), PARAMETER :: io_ace = 11, &   ! input file (used to fixed with name "model_setup_BTE.in")
                               io_debug = 12, & ! file for debug printout
                               io_user_input = 13, & ! file with input info for simulation ("user_polar.in")
                               io_residual= 14, & ! file with  count and residual values
                               io_output = 15, &     ! file for temperature contour
                               io_flux = 16  ! File for heat flux output
  
  CHARACTER(len=32) :: nProcsName, rankName
  LOGICAL, PARAMETER :: print_tnot = .true. 
  LOGICAL, PARAMETER :: print_tempInterface = .false.
  LOGICAL, PARAMETER :: print_fluxInterface = .false.
  LOGICAL, PARAMETER :: print_tbcInterface = .false.
  LOGICAL, PARAMETER :: print_interiorFlux = .false.
  INTEGER, PARAMETER :: io_tnot = 94, io_tbcInterface = 96, io_interiorFlux = 31, io_tnot_vtk_series = 95, io_tnot_vtk=93

  !...hard-code to choose cell type for plotting the mesh in vtk format
  !...currently support: 10 (VTK_TETRA), 12 (VTK_HEXAHEDRON), 13 (VTK_EDGE), 9 (VTK_QUAD)
  INTEGER, PARAMETER :: vtk_cellType = 10
  
  ! Debug level indicates granularity of printing
  ! 0 => no printing
  ! 5 => everything printed
  INTEGER(int_p), PARAMETER :: debug_level = 0, &   !...0 is no printing to debug.out
                              max_iter = 3000, &   !...maximum number of outer iterations
				                      min_iter = 1, & !5, &       !...minimum number of outer iterations	
                              iter_gmres = 100, &   !...maximum number of GMRES iterations
                              freq_time_output = 1  !...output frequency
  !! la and ta phonons                             
  INTEGER(int_p), PARAMETER :: LA_1 = 1, TA_1 = 2 !TA_1= 3   ,&      !silicon
                               !LA_2 = 2, TA_2 =4       !! germanium
                               
  INTEGER(int_p), PARAMETER :: SILICON = 1, GERMENIUM = 2, ALUMINUM = 3
  
  REAL(real_p), PARAMETER :: relax_temp = 1.0d0 !, & ! Relaxation factor for Temperature
                             !tol_inner = 1.0D-2 ! , & ! GMRES tolerance
                             !tol_outer = 1.0D-9  ! Outer loop tolerance
  REAL(real_p) :: tol_outer

  REAL(4) :: time_start,time_end,total_time,ta(2) 
  REAL(real_p) :: residual,dt,idt 
  REAL(real_p) :: residual_sum !...sum of residual across ranks
  REAL(real_p) :: dtheta,dphi,hotwall,coldwall
  REAL(real_p) :: initial_temp, dos
  REAL(real_p) :: qdblprime, laserpower, spotsize, qnot, layer
  REAL(real_p) :: qnot_g 
  !REAL(real_p):: delta_la,delta_ta
  
  !...REAL(real_p), Allocatable, Dimension(:) :: alu
  REAL(real_p), Allocatable, Dimension(:) :: rmu,rxi,ret,inrmu,inrxi,inret
  REAL(real_p), Allocatable, Dimension(:) :: Omega
  !  REAL(real_p), Allocatable, Dimension(:) :: gpvel
  REAL(real_p), Allocatable, Dimension(:) :: centfreq
  !...REAL(real_p), Allocatable  ::  Am(:)
  !...REAL(real_p), Allocatable  ::  R(:), soln(:)
  REAL(real_p), Dimension(:), Pointer :: tnot, tprime, tnotv, tbprime, tnotone, &
    temp_bcone, ttop, tbot, ttopone  ! flux
  REAL(real_p), Dimension(:), Pointer :: s
  !...REAL(real_p), Dimension(:), Pointer :: resid
  !...REAL(real_p), Dimension(:), Pointer :: aps !, scs_tran      !!!!!!scs_tran added 28th aug
  !...REAL(real_p), Dimension(:,:) , Pointer :: anbs
  !...REAL(real_p), Allocatable, Dimension(:,:) :: scs
  !...REAL(real_p), Dimension(:,:) , Pointer ::  scs_tran
  REAL(real_p), Dimension(20):: xi,wi

  !REAL(real_p), DIMENSION(:,:), POINTER :: sdotn_array, insdotn_array !...Feb 22, 2022: pre-computed sdotn and insdotn
  
  REAL(real_p), Dimension(:,:) , Pointer ::  wavenumber
  !...REAL(real_p), Allocatable, Dimension(:,:) :: vv
  
  REAL(real_p), Dimension(:,:,:), Pointer :: gc, gcone, flux, jfacepos, jfaceneg , inot_c, &
    flux_s, jfacepos_s, jfaceneg_s !!! *_s for send buffer.
  REAL(real_p), Dimension(:,:,:,:), Pointer :: intensity, intensityone

  !! interior flux
  !...Feb 4, 2022: no more computation of interior flux
  !REAL(real_p), DIMENSION(:,:) , Pointer :: flux_crossfar, flux_tranfar, flux_crossnear, flux_trannear
  !ßREAL(real_p), DIMENSION(:,:) , Pointer :: flux_crossfar_s, flux_tranfar_s, flux_crossnear_s, flux_trannear_s

  REAL(real_p), ALLOCATABLE, DIMENSION(:) :: delta_band, band_low
  REAL(real_p), ALLOCATABLE, DIMENSION(:,:) :: gpvel
  REAL(real_p), ALLOCATABLE, DIMENSION(:) :: wdisper

  REAL(real_p), ALLOCATABLE, DIMENSION(:) :: qgen

  REAL(real_p), ALLOCATABLE, DIMENSION(:) :: k_sil
  
  !!!!! for variable dt
  INTEGER (int_p) ::  timestep_period, n           !!  no of pulse!
  !! for dt calculation
  INTEGER(int_p) :: ntimep, npulse !! ntimep = no of timestep in pluse, npulse== no of pulse
  REAL(real_p) :: dt_1  !! first dt of the pulse
  REAL(real_p) :: time_pulse
  REAL(real_p) :: sep !! will be calculated

  !!!!! for harmonic q

  REAL(real_p) :: mod_freq
  REAL(real_p) :: rel
  REAL(real_p) :: probe_size !...for ALSI boundary condition
  REAL(real_p) :: weighted_ttop_local_sum, bfArea_local_sum, weighted_ttop_global_sum, bfArea_global_sum !...for ALSI boundary condition
End Module variables
