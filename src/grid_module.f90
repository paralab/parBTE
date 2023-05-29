!***********************************************************************
! Parameters and data used ! 3d grid

MODULE GRID

  Use PRECISIONS
  Implicit None
  SAVE
  ! Declare global variables
  ! USER GLOBAL VARIABLE DECLARATION BEGIN
  ! number of cells, faces etc.
  INTEGER(int_p) :: ncells, nfaces, nnodes, nbcfaces
  INTEGER :: nf_max, nbc_patches       
  INTEGER(int_p) :: alsi_count, axi_check !...alsi_count is used in aluminiumtemp_ltc_exp.f90
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: face_alsi !...is used in aluminiumtemp_ltc_exp.f90
  INTEGER(int_p), PARAMETER :: ADIA = 4, ISOT = 3, CONJ = 8, HFLUX = 5, SYMM = 2
  INTEGER(int_p), PARAMETER :: ALSI = 9             !!! taking the al/si boundary as thinwall from model file.
  ! arrays relating to boundary faces.
  ! bface stores bface index
  ! f_to_bf is mapping from global face index to bface
  ! bf_to_f is mapping from bface to global face index
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: f_to_bf, bf_to_f, bf_to_c, &
                                              nfcell, bctype, bcname, nvcell , cell_type, material_type
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: bnode, bface !...Han: I think bnode is never used ?
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: nvface, no_f_patch 
  INTEGER(int_p), DIMENSION(:,:), ALLOCATABLE :: normdir 
  INTEGER(real_p), DIMENSION(:,:), ALLOCATABLE :: lcf, lcv, lfv, lfc, lcc 

  ! cell center coordinates, face center coordinates, vertex coordinates, face area
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: xc, yc, zc, xf, yf, zf, xv, yv, zv, areaf, &
                                             vecfx, vecfy, vecfz, volcell
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: disn, phibc, wcf 
  REAL(real_p), DIMENSION(:,:), ALLOCATABLE :: wcv, wbfv  
  Real(real_p), DIMENSION(:), ALLOCATABLE :: temp_bc
  
  LOGICAL :: threed = .true.

  !...Han added, for cell-based parallel
  !...local cell id to global cell id, lc_to_gc(1:ncells)--> owned cells; lc_to_gc(ncells+1:nCells_total)--> ghost cells
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: lc_to_gc
  !...Aug 06: map from local cell id to global cell id sorted by partition, i.e. cells=[1,...,n0] have partion 0, cells=[n0+1,...,n1] have partition 1
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: lc_to_pgc
  !...number of ghost cells, nCells_total = ncells + nGhostCells
  INTEGER(int_p) :: nGhostCells, nCells_total
  
END MODULE GRID