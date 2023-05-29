! ***************************************************************!
! This subroutine reads in all data written out from 
! CFD-ACE+ (using ace_output.f) into the file model_setup_BTE.in
! It can handle any mesh both in 2D and 3D
! The subroutines processes all the information and finally
! fills in the relevant arrays in this code.
!
! ***************************************************************!
! \brief read in mesh file and construct local mesh
! \authors Han D. Tran, Siddharth Saurav, P. Sadayappan, Sandip Mazumder, Hari Sundar
! \date 2021-2023
! \details
! subroutine to compute average ttop for boundary faces within the probe area, applied for ALSI boundary condition
! version 03: lfc, lcc, bf_to_c gives local cell id,
! bface(i) = 2 if face i is interface with neighbor
! required: global cell ids are contiguous:
! ghost cells of global ids < my global-id-range are owned by ranks < myRank
! ghost cells of global ids > my global-id-range are owned by ranks > myRank
SUBROUTINE grid_and_boundary(modelFile)

  USE precisions
  USE grid
  USE constants, ONLY: zero, one
  USE variables, ONLY: io_ace, io_debug, debug_level, rank, rank_cart, rank_cell, np_max, num_material

  IMPLICIT NONE

  !...Han added model filename get from runtime
  CHARACTER(256), INTENT(IN) :: modelFile

  INTEGER(int_p) :: icell, iface, inode, ibf, i, j, k, cell_1, cell_2, node_1, &
      node_2, ifc, ic, ibc, throw_away1, throw_away2
  INTEGER(int_p) :: nv_max1, nv_max2, ierr, flag2
  
  REAL(real_p) :: xtemp1, xtemp2, disntemp
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: sumweight, sumwtf
  
  LOGICAL :: exists = .true.
  LOGICAL :: flag1 = .false.

  !...variables holding global grid (same for all ranks), just local to this subroutine
  INTEGER(int_p) :: ncells_GB, nfaces_GB, nnodes_GB, nbcfaces_GB
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: nfcell_GB, nvcell_GB, owner_GB
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: gc_to_pgc !...Aug 6: map from original global cell id to petsc global id that satisfies the contiguours labeling condition
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: pgc_to_gc !...Aug 23: map from petsc global cell id to original global cell id
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: xc_GB, yc_GB, zc_GB, volcell_GB
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: xf_GB, yf_GB, zf_GB
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: vecfx_GB, vecfy_GB, vecfz_GB, areaf_GB
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: xv_GB, yv_GB, zv_GB

  INTEGER(int_p), DIMENSION(:,:), ALLOCATABLE :: lcf_GB, lcv_GB, lfc_GB, lfv_GB, lcc_GB
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: nvface_GB
  INTEGER(int_p) :: nbc_patches_GB
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: no_f_patch_GB !...no_f_patch declared in GRID is no longer used, could be deleted
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: phibc_GB !...phibc is only used in this subroutine, could be deallocated at the end of this subroutine
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: bface_GB, f_to_bf_GB
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: bf_to_f_GB, bf_to_c_GB, bctype_GB
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: cell_type_GB, material_type_GB

  INTEGER(int_p) nfaces_upper_bound, nnodes_upper_bound
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: lf_to_gf, lv_to_gv, lbf_to_gbf
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: gc_to_lc, gf_to_lf, gv_to_lv,  gbf_to_lbf
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: oc_to_gc ! map from owned cells to global id
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: ghost_cells_gc !list of (repeated) global id of ghost cells
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: ghost_xc, ghost_yc, ghost_zc !for computing disn(), wcf()
  INTEGER(int_p), DIMENSION(:), ALLOCATABLE :: ghost_cells_pgc !list of ghost cells in terms of petsc global ids

  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug, '(a)') "Starting grid_and_boundary" 
  
  INQUIRE(file = modelFile, EXIST = exists)
  IF (.not. exists) THEN
    WRITE(*,'(a)') "Model file could not be found!"
    STOP
  ENDIF

  OPEN (UNIT=io_ace, FILE=modelFile, status = "old")

  ! Reading in global grid size data
  READ(io_ace,*) ncells_GB, nfaces_GB, nnodes_GB, nbcfaces_GB
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug, '(a)') &
    "Finished reading global grid size data"

  ! Read in other options
  !READ(io_ace,*) threed,flag1
  READ(io_ace,*) throw_away1, throw_away2
  IF (throw_away1 .eq. 0) THEN
    threed = .false.
  ELSE
    threed = .true.
  ENDIF

  IF (throw_away2 .eq. 0) THEN
    flag1 = .false.
  ELSE
    flag1 = .true.              ! SOLID PRESENT
  ENDIF

  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a,2(1x,i0))') &
    "Finished reading other options", threed, flag1

  ! Reading in number of faces and vertices, and owner for each cell
  ALLOCATE(nfcell_GB(ncells_GB), stat=ierr)
  ALLOCATE(nvcell_GB(ncells_GB), stat=ierr)
  ALLOCATE(owner_GB(ncells_GB), stat=ierr)
  DO i = 1,ncells_GB
    READ(io_ace,*) nfcell_GB(i), nvcell_GB(i), owner_GB(i)
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
          "Finished reading number of faces, vertices, and owner for each cell"
  
  ! Reading in geometric data
  IF (threed) THEN
    !...only for owned cell: coordinates of cell center, cell volume
    ALLOCATE (xc_GB(ncells_GB), stat=ierr)
    ALLOCATE (yc_GB(ncells_GB), stat=ierr)
    ALLOCATE (zc_GB(ncells_GB), stat=ierr)
    ALLOCATE (volcell_GB(ncells_GB), stat=ierr)

    DO i = 1,ncells_GB
      READ(io_ace,*) xc_GB(i), yc_GB(i), zc_GB(i), volcell_GB(i)
    END DO

    ALLOCATE (xf_GB(nfaces_GB), stat=ierr)
    ALLOCATE (yf_GB(nfaces_GB), stat=ierr)
    ALLOCATE (zf_GB(nfaces_GB), stat=ierr)
    ALLOCATE (vecfx_GB(nfaces_GB), stat=ierr)
    ALLOCATE (vecfy_GB(nfaces_GB), stat=ierr)
    ALLOCATE (vecfz_GB(nfaces_GB), stat=ierr)
    ALLOCATE (areaf_GB(nfaces_GB), stat=ierr)

    DO i = 1,nfaces_GB
      READ(io_ace,*) xf_GB(i), yf_GB(i), zf_GB(i), &
        areaf_GB(i),vecfx_GB(i),vecfy_GB(i),vecfz_GB(i)
    END DO
    
    ALLOCATE (xv_GB(nnodes_GB), stat=ierr)
    ALLOCATE (yv_GB(nnodes_GB), stat=ierr)
    ALLOCATE (zv_GB(nnodes_GB), stat=ierr)

    DO i = 1,nnodes_GB
      READ(io_ace,*) xv_GB(i), yv_GB(i), zv_GB(i)
    END DO
  ELSE
    ALLOCATE (xc_GB(ncells_GB), stat=ierr)
    ALLOCATE (yc_GB(ncells_GB), stat=ierr)
    ALLOCATE (volcell_GB(ncells_GB), stat=ierr)

    DO i = 1,ncells_GB
      READ(io_ace,*) xc_GB(i), yc_GB(i), volcell_GB(i)
    END DO

    ALLOCATE (xf_GB(nfaces_GB), stat=ierr)
    ALLOCATE (yf_GB(nfaces_GB), stat=ierr)
    ALLOCATE (vecfx_GB(nfaces_GB), stat=ierr)
    ALLOCATE (vecfy_GB(nfaces_GB), stat=ierr)
    ALLOCATE (areaf_GB(nfaces_GB), stat=ierr)

    DO i = 1,nfaces_GB
      READ(io_ace,*) xf_GB(i), yf_GB(i), areaf_GB(i), vecfx_GB(i), vecfy_GB(i)
    END DO

    ALLOCATE (xv_GB(nnodes_GB), stat=ierr)
    ALLOCATE (yv_GB(nnodes_GB), stat=ierr)

    DO i = 1,nnodes_GB
      READ(io_ace,*) xv_GB(i), yv_GB(i)
    END DO
  ENDIF
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading geometric data"

  ! Reading cell-to-face connectivity
  !...determine max n faces/cell (only search on cells owned by rank)
  nf_max = -1
  DO i = 1,ncells_GB
    nf_max = MAX(nf_max, nfcell_GB(i))
  END DO

  ALLOCATE(lcf_GB(ncells_GB, nf_max), stat=ierr)

  DO i = 1,ncells_GB
    READ(io_ace,*) (lcf_GB(i,j), j = 1,nfcell_GB(i))
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading cell to face connectivity"

  ! Reading cell-to-vertex connectivity
  nv_max1 = -1
  DO i = 1,ncells_GB
    nv_max1 = MAX(nv_max1, nvcell_GB(i))
  END DO

  !...cell-to-vertex map
  ALLOCATE(lcv_GB(ncells_GB, nv_max1), stat=ierr)
  
  DO i = 1,ncells_GB
    READ(io_ace,*) (lcv_GB(i,j), j = 1,nvcell_GB(i))
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading cell to vertex connectivity"

  ! Reading face-to-cell connectivity
  ALLOCATE(lfc_GB(nfaces_GB, 2), stat=ierr)
  DO i = 1,nfaces_GB
    READ(io_ace,*) lfc_GB(i,1), lfc_GB(i,2)
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading face to cell connectivity"

  ! Reading face-to-vertex connectivity
  ALLOCATE(nvface_GB(nfaces_GB), stat=ierr)
  nv_max2 = -1
  DO i = 1,nfaces_GB
    READ(io_ace,*) nvface_GB(i)
    nv_max2 = MAX(nv_max2, nvface_GB(i))
  END DO

  ALLOCATE(lfv_GB(nfaces_GB, nv_max2), stat=ierr)
  
  DO i = 1,nfaces_GB
    READ(io_ace,*) (lfv_GB(i,j), j = 1,nvface_GB(i))
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading face to vertex connectivity"

  ! Reading cell-to-cell connectivity
  ALLOCATE(lcc_GB(ncells_GB, nf_max), stat=ierr)
  DO i = 1,ncells_GB
    READ(io_ace,*) (lcc_GB(i,j), j=1,nfcell_GB(i))
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading cell type and cell to cell connectivity"


  ! Reading in boundary data
  READ(io_ace,*) nbc_patches_GB

  ALLOCATE(no_f_patch_GB(nbc_patches_GB), stat=ierr)  ! No of faces per patch
  ALLOCATE(phibc_GB(nbcfaces_GB), stat=ierr)
      
  ibf = 0 
  DO j = 1, nbc_patches_GB
    READ(io_ace,*) no_f_patch_GB(j)
    DO i = 1, no_f_patch_GB(j)
      ! Reading values of all variables at boundaries
      ibf = ibf + 1
      READ(io_ace,*) phibc_GB(ibf)
    ENDDO
  ENDDO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading boundary data"

  ! Reading Boundary face flag and global face to boundary face conn
  ALLOCATE(bface_GB(nfaces_GB), stat=ierr)
  ALLOCATE(f_to_bf_GB(nfaces_GB), stat=ierr)

  DO i = 1,nfaces_GB
    READ(io_ace,*) bface_GB(i), f_to_bf_GB(i)
  END DO

  ! Reading Boundary type and boundary face to global face and cell
  ALLOCATE(bf_to_f_GB(nbcfaces_GB), stat=ierr)
  ALLOCATE(bf_to_c_GB(nbcfaces_GB), stat=ierr)
  ALLOCATE(bctype_GB(nbcfaces_GB), stat=ierr)
  DO i = 1,nbcfaces_GB
    READ(io_ace,*) bctype_GB(i), bf_to_f_GB(i), bf_to_c_GB(i)
  END DO

  ! Reading cell material data   added on 16th april, not readin material type
  ! IF (flag1) THEN
    ALLOCATE(cell_type_GB(ncells_GB), stat=ierr)
    DO i = 1,ncells_GB
      READ(io_ace,*)  cell_type_GB(i)
    ENDDO

    ALLOCATE(material_type_GB(ncells_GB), stat=ierr)
    DO i = 1,ncells_GB
      READ(io_ace,*)  material_type_GB(i)
    ENDDO
  ! ENDIF

  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished reading boundary connectivity"

  !READ(io_ace,*) axi_check
  num_material = 1
  np_max = num_material * 2

  !...extract local data from global data ============================================================
  
  !...Aug 06: based on partition, determine the global ID that satisfy the contiguous labeling condition
  ALLOCATE(gc_to_pgc(ncells_GB), stat=ierr)
  DO i = 1,ncells_GB
    k = 1
    DO j = 1,ncells_GB
      IF (owner_GB(j) < owner_GB(i)) k = k+1
      IF ((owner_GB(j).eq.owner_GB(i)).and.(j.lt.i)) k = k+1
    ENDDO
    gc_to_pgc(i) = k
  ENDDO

  !...Aug 23: map from petsc global cell id to original global cell id (for use in ghost cells...)
  ALLOCATE(pgc_to_gc(ncells_GB), stat=ierr)
  DO i = 1,ncells_GB
    pgc_to_gc(gc_to_pgc(i)) = i
  ENDDO

  !...map from global cell/face/vertex id to local id (this could be too expensive for memory!)
  ALLOCATE(gc_to_lc(ncells_GB), stat=ierr)
  ALLOCATE(gf_to_lf(nfaces_GB), stat=ierr)
  ALLOCATE(gv_to_lv(nnodes_GB), stat=ierr)
  ALLOCATE(gbf_to_lbf(nbcfaces_GB), stat=ierr)
  gc_to_lc = -1
  gf_to_lf = -1
  gv_to_lv = -1
  gbf_to_lbf = -1

  !...determine number of owned cells
  !...also determine upper bounds (i.e. including redundants) for local number of faces & vertices
  ncells = 0
  nfaces_upper_bound = 0
  nnodes_upper_bound = 0
  DO i = 1,ncells_GB
    IF (owner_GB(i).eq.rank_cell) THEN
      ncells = ncells + 1
      nfaces_upper_bound = nfaces_upper_bound + nfcell_GB(i) !...TODO: this could be too large number causing insufficient memory
      nnodes_upper_bound = nnodes_upper_bound + nvcell_GB(i) !...TODO: this could be too large number causing insufficient memory
    ENDIF
  ENDDO

  !...with ncells determined, allocate: xc, yc, zc, nfcell, nvcell, volcell, lcf, lcv, lcc
  ALLOCATE(nfcell(ncells), stat=ierr)
  ALLOCATE(nvcell(ncells), stat=ierr)
  ALLOCATE(xc(ncells), stat=ierr)
  ALLOCATE(yc(ncells), stat=ierr)
  IF (threed) THEN
    ALLOCATE(zc(ncells), stat=ierr)
  ENDIF
  ALLOCATE(volcell(ncells), stat=ierr)
  ALLOCATE(lcf(ncells, nf_max), stat=ierr)
  ALLOCATE(lcv(ncells, nv_max1), stat=ierr) !nv_max1 is max vertices per cell (nv_max2 is max vertices per face)
  !ALLOCATE(lcc(ncells, nf_max), stat=ierr)
  ALLOCATE(lcc(nf_max, ncells), stat=ierr) !...2021.09.13

  !...allocate local-to-global maps of faces & vertices
  ALLOCATE(lf_to_gf(nfaces_upper_bound), stat=ierr)
  ALLOCATE(lv_to_gv(nnodes_upper_bound), stat=ierr)

  ALLOCATE(oc_to_gc(ncells), stat=ierr) !...map from owned cell to global cell id

  !...determine number of local faces & vertices
  icell = 0
  iface = 0
  inode = 0
  DO i = 1,ncells_GB
    IF (owner_GB(i).eq.rank_cell) THEN
      !...local cell id
      icell = icell + 1
      !...map of local cell id to global cell id
      oc_to_gc(icell) = i
      !...cell center of local cell
      xc(icell) = xc_GB(i)
      yc(icell) = yc_GB(i)
      IF (threed) THEN
        zc(icell) = zc_GB(i)
      ENDIF
      !...cell volume of local cell
      volcell(icell) = volcell_GB(i)
      !...number of faces & vertices of local cell
      nfcell(icell) = nfcell_GB(i)
      nvcell(icell) = nvcell_GB(i)

      !...push all faces composing owned cells to the (pre-sorted and redundant) map
      DO j = 1,nfcell_GB(i)
        iface = iface + 1
        lf_to_gf(iface) = lcf_GB(i,j)
      ENDDO

      !...push all vertices composing owned cells to local-to-global maps (pre-sorted and contain redundants)
      DO j = 1,nvcell_GB(i)
        inode = inode + 1
        lv_to_gv(inode) = lcv_GB(i,j)
      ENDDO

    ENDIF ! owner_GB(i) .eq. rank_cell
  ENDDO ! i = 1,ncells_GB

  !...build global-to-local cell map (for owned cells, later will include ghost cells)
  DO icell = 1,ncells
    gc_to_lc(oc_to_gc(icell)) = icell
  ENDDO

  !...sorting map from local-to-global face/vertex id
  !...this prepares for naming local face/vertex id as in the order of global face/vertex that composing celles owned by me
  CALL quicksort_int(lf_to_gf, 1, nfaces_upper_bound)
  CALL quicksort_int(lv_to_gv, 1, nnodes_upper_bound)

  !...delete redundant face ids in lf_to_gf, then determine (exact) number of faces
  k = 1
  DO i = 2, nfaces_upper_bound
    IF (lf_to_gf(i) .ne. lf_to_gf(k)) THEN
      k = k + 1
      lf_to_gf(k) = lf_to_gf(i)
    ENDIF
  ENDDO
  !...update number of local cells
  nfaces = k

  !...delete redundant vertex ids in lv_to_gv, then determine (exact) number of vertices
  k = 1
  DO i = 2, nnodes_upper_bound
    IF (lv_to_gv(i) .ne. lv_to_gv(k)) THEN
      k = k + 1
      lv_to_gv(k) = lv_to_gv(i)
    ENDIF
  ENDDO
  !...update number of local vertices
  nnodes = k

  !...build global-to-local maps for face & vertex
  DO iface = 1,nfaces
    gf_to_lf(lf_to_gf(iface)) = iface
  ENDDO
  DO inode = 1,nnodes
    gv_to_lv(lv_to_gv(inode)) = inode
  ENDDO

  !...compute geometric data of faces & vertices: xf, yf, zf, vecfx, vecfy, vecfz, areaf, xv, yv, zv
  ALLOCATE(xf(nfaces), stat=ierr)
  ALLOCATE(yf(nfaces), stat=ierr)
  ALLOCATE(vecfx(nfaces), stat=ierr)
  ALLOCATE(vecfy(nfaces), stat=ierr)
  ALLOCATE(areaf(nfaces), stat=ierr)
  ALLOCATE(xv(nnodes), stat=ierr)
  ALLOCATE(yv(nnodes), stat=ierr)
  IF (threed) THEN
    ALLOCATE(zf(nfaces), stat=ierr)
    ALLOCATE(vecfz(nfaces), stat=ierr)
    ALLOCATE(zv(nnodes), stat=ierr)
  ENDIF
  
  DO iface = 1,nfaces
    xf(iface) = xf_GB(lf_to_gf(iface))
    yf(iface) = yf_GB(lf_to_gf(iface))
    vecfx(iface) = vecfx_GB(lf_to_gf(iface))
    vecfy(iface) = vecfy_GB(lf_to_gf(iface))
    areaf(iface) = areaf_GB(lf_to_gf(iface))
    IF (threed) THEN
      zf(iface) = zf_GB(lf_to_gf(iface))
      vecfz(iface) = vecfz_GB(lf_to_gf(iface))
    ENDIF
  ENDDO

  DO inode = 1,nnodes
    xv(inode) = xv_GB(lv_to_gv(inode))
    yv(inode) = yv_GB(lv_to_gv(inode))
    IF (threed) THEN
      zv(inode) = zv_GB(lv_to_gv(inode))
    ENDIF
  ENDDO

  !...determine lcf and lcv in terms of local face id and local vertex id
  DO icell = 1,ncells
    !...lcf
    DO j = 1,nfcell(icell)
      lcf(icell, j) = gf_to_lf(lcf_GB(oc_to_gc(icell), j))
    ENDDO
    !...lcv
    DO j = 1,nvcell(icell)
      lcv(icell, j) = gv_to_lv(lcv_GB(oc_to_gc(icell), j))
    ENDDO
  ENDDO !icell = 1,ncells

  !...allocate variables holding boundary faces
  ALLOCATE(bface(nfaces), stat=ierr)
  ALLOCATE(f_to_bf(nfaces), stat=ierr)
  f_to_bf = 0

  !...determine local number of boundary faces
  nbcfaces = 0
  DO iface = 1,nfaces
    bface(iface) = bface_GB(lf_to_gf(iface))      !flag if iface is boundary or not
    IF (bface(iface) .eq. 1) THEN
      nbcfaces = nbcfaces + 1
    ENDIF
    f_to_bf(iface) = f_to_bf_GB(lf_to_gf(iface))  !this is still the global bdr face id
  ENDDO

  ALLOCATE(lbf_to_gbf(nbcfaces), stat=ierr)

  ALLOCATE(bctype(nbcfaces), stat=ierr)
  ALLOCATE(bf_to_f(nbcfaces), stat=ierr)
  ALLOCATE(bf_to_c(nbcfaces), stat=ierr)

  ALLOCATE(phibc(nbcfaces), stat=ierr)

  !...build map from local bf to global bf
  k = 0
  DO iface = 1,nfaces
    IF (bface(iface).eq.1) THEN
      k = k + 1
      lbf_to_gbf(k) = f_to_bf(iface) !remind: f_to_bf() is still global bf id
    ENDIF
  ENDDO
  !...sort so that local bf is in order of global bf
  CALL quicksort_int(lbf_to_gbf, 1, nbcfaces)

  !...build map from global bf to local bf
  DO iface = 1,nbcfaces
    gbf_to_lbf(lbf_to_gbf(iface)) = iface
  ENDDO

  !...update map from local f to local bf
  DO iface = 1,nfaces
    IF (bface(iface).eq.1) THEN
      f_to_bf(iface) = gbf_to_lbf(f_to_bf(iface))
    ENDIF
  ENDDO
  
  !...update corresponding variables
  DO iface = 1,nbcfaces
    bctype(iface) = bctype_GB(lbf_to_gbf(iface))
    bf_to_f(iface) = gf_to_lf(bf_to_f_GB(lbf_to_gbf(iface)))
    bf_to_c(iface) = gc_to_lc(bf_to_c_GB(lbf_to_gbf(iface)))
    phibc(iface) = phibc_GB(lbf_to_gbf(iface))
  ENDDO

  !...update cell type and material type
  ALLOCATE(cell_type(ncells), stat=ierr)
  ALLOCATE(material_type(ncells), stat=ierr)
  DO i = 1,ncells
    cell_type(i) = cell_type_GB(oc_to_gc(i))
    material_type(i) = material_type_GB(oc_to_gc(i))
  ENDDO

  !...determine ghost cells
  !...so far we do not know exactly how many ghost cells, thus allocate with max ghost cells = nfaces - nbcfaces
  ALLOCATE(ghost_cells_gc(nfaces - nbcfaces), stat=ierr)

  !...loop over faces, if an associated cell is not a local cell then add it to the list of ghost cells
  k = 0
  DO iface = 1,nfaces
    !...boundary face cannot associated with ghost cell
    IF (bface(iface).eq.1) THEN
      CYCLE
    ENDIF
    !...check if side-1 cell is a ghost cell
    IF (gc_to_lc(lfc_GB(lf_to_gf(iface),1)) .eq. -1) THEN
      k = k + 1
      ghost_cells_gc(k) = lfc_GB(lf_to_gf(iface),1)
    ELSE
      !...if side-1 cell is an owned cell, continue to check if side-2 cell is a ghost
      IF (gc_to_lc(lfc_GB(lf_to_gf(iface),2)) .eq. -1) THEN
        k = k + 1
        ghost_cells_gc(k) = lfc_GB(lf_to_gf(iface),2)
      ENDIF
    ENDIF
  ENDDO

  !...sorting (only sort the first k components because the remaining is junk)
  !...note: this sorting is on initial global id (not "pgc" which is used in both matrix-assembled and assembly-free)
  IF (k .ge. 1) THEN
    CALL quicksort_int(ghost_cells_gc, 1, k)
  ENDIF

  !...remove redundants
  !...bug 2021.08.03, only do this if k >= 1 (at least there is one ghost cells)
  IF (k .ge. 1) THEN
    j = 1
    DO i = 2, k
      IF (ghost_cells_gc(i) .ne. ghost_cells_gc(j)) THEN
        j = j + 1
        ghost_cells_gc(j) = ghost_cells_gc(i)
      ENDIF
    ENDDO
    !...update number of ghost cells
    nGhostCells = j
  ENDIF

  !...total number of cells including ghosts
  ncells_total = ncells + nGhostCells

  !... 2021.08.23: sorting ghost cells in terms of "pgc"
  ALLOCATE(ghost_cells_pgc(nGhostCells), stat=ierr)
  DO i = 1,nGhostCells
    ghost_cells_pgc(i) = gc_to_pgc(ghost_cells_gc(i))
  ENDDO
  
  IF (nGhostCells .ge. 1) THEN
    CALL quicksort_int(ghost_cells_pgc, 1, nGhostCells)
  ENDIF

  !...after sorting ghost pgc, we need to update ghost_cells_gc
  !...eg.: 2D_5by5_p5b.in, for P1:
  !...ghost_cells_gc = [3,8,12,16] corresponds to (before sorting): ghost_cells_pgc = [11,14,1,2]
  !...then, after sorting: ghost_cells_pgc = [1,2,11,14]
  !...thus, we need to update as: ghost_cells_gc = [pgc_to_gc(1),pgc_to_gc(2),pgc_to_gc(11),pgc_to_gc(14)] = [12,16,3,8]
  DO i = 1,nGhostCells
    ghost_cells_gc(i) = pgc_to_gc(ghost_cells_pgc(i))
  ENDDO
  
  !...map from local cells (including ghost celss) to global cells
  ALLOCATE(lc_to_gc(ncells_total), stat=ierr) !this belongs to GRID (declared in grid_module.f90)
  DO icell = 1,ncells
    lc_to_gc(icell) = oc_to_gc(icell)
  ENDDO
  !...eg, 2D_5by5_p5b.in, for P1, after this we have
  !...lc_to_gc = [1,2,6,7,11,12,16,3,8] in which 12,16,3,8 are ghost cells
  DO icell = 1,nGhostCells
    !...2021.08.23 change local id for ghost cells according to pgc
    !lc_to_gc(ncells + icell) = ghost_cells_gc(icell)
    lc_to_gc(ncells + icell) = pgc_to_gc(ghost_cells_pgc(icell))
  ENDDO
  
  !...Aug 06, add lc_to_pgc(i) from global cell id to the sorted cell id that satisfy the contiguous labeling condition
  ALLOCATE(lc_to_pgc(ncells_total), stat=ierr)
  DO icell = 1,ncells+nghostcells
    lc_to_pgc(icell) = gc_to_pgc(lc_to_gc(icell))
  ENDDO
  
  !...update gc_to_lc including ghost cells
  !...eg, 2D_5by5_p5b.in, for P1, after this we have:
  !...gc_to_lc(12)=6, gc_to_lc(16)=7 (Note: before when we sort ghost cells in term of gc, this was gc_to_lc(12)=8 and gc_to_lc(16)=9)
  !...basically, we want to put ghost cells in the sorted order of "pgc" instead of gc because this gives correct communication used in the aMat approach
  DO icell = 1,nghostcells
    !...2021.08.23 change local id for ghost cells according to pgc
    !gc_to_lc(ghost_cells_gc(icell)) = ncells + icell
    gc_to_lc(pgc_to_gc(ghost_cells_pgc(icell))) = ncells + icell
  ENDDO

  !...lfc: two (LOCAL) cell ids associated with (local) face
  ALLOCATE(lfc(nfaces,2), stat=ierr)
  DO iface = 1,nfaces
    lfc(iface,1) = gc_to_lc(lfc_GB(lf_to_gf(iface),1))
    lfc(iface,2) = gc_to_lc(lfc_GB(lf_to_gf(iface),2))
  ENDDO
  
  !...lfv: local vertex ids composing (local) face
  ALLOCATE(nvface(nfaces), stat=ierr)
  ALLOCATE(lfv(nfaces, nv_max2), stat=ierr)
  DO iface = 1,nfaces
    nvface(iface) = nvface_GB(lf_to_gf(iface))
    DO j = 1,nvface(iface)
      lfv(iface, j) = gv_to_lv(lfv_GB(lf_to_gf(iface), j))
    ENDDO
  ENDDO
  
  !...lcc: (LOCAL) cell ids associated with faces composing (local) cell
  DO icell = 1,ncells
    DO j = 1,nfcell(icell)
      !lcc(icell, j) = gc_to_lc(lcc_GB(lc_to_gc(icell), j))
      lcc(j,icell) = gc_to_lc(lcc_GB(lc_to_gc(icell), j))
    ENDDO
  ENDDO

  !...cell center of ghost cells, for computing disn and wcf
  ALLOCATE(ghost_xc(nghostcells), stat=ierr)
  ALLOCATE(ghost_yc(nghostcells), stat=ierr)
  IF (threed) THEN
    ALLOCATE(ghost_zc(nghostcells), stat=ierr)
  ENDIF
  DO i = 1,nghostcells
    ghost_xc(i) = xc_GB(ghost_cells_gc(i))
    ghost_yc(i) = yc_GB(ghost_cells_gc(i))
    IF (threed) THEN
      ghost_zc(i) = zc_GB(ghost_cells_gc(i))
    ENDIF
  ENDDO
  !==========================================================================================
  ! Finished all Reading. Now time to process data
  ! Process Distance
  !... I do not see disn is used anywhere else in the code ?
  ALLOCATE (disn(nfaces), stat=ierr)
  IF(threed)THEN
    DO iface=1,nfaces
        cell_1 =  lfc(iface,1)
        cell_2 =  lfc(iface,2)
        IF(cell_1 == cell_2) THEN  ! External faces
          !...boundary face cannot be associated with a ghost cell
          disn(iface) = ABS((xf(iface) - xc(cell_1))*vecfx(iface) + &
                            (yf(iface) - yc(cell_1))*vecfy(iface) + &
                            (zf(iface) - zc(cell_1))*vecfz(iface))
        ELSE  ! Interior face
          IF (cell_1 > ncells) THEN
            !...cell_1 is ghost cell, then cell_2 cannot be ghost
            disn(iface) = ABS((xc(cell_2) - ghost_xc(cell_1 - ncells))*vecfx(iface) + &
                            (yc(cell_2) - ghost_yc(cell_1 - ncells))*vecfy(iface) + &
                            (zc(cell_2) - ghost_zc(cell_1 - ncells))*vecfz(iface))
          ELSEIF (cell_2 > ncells) THEN
            !...cell_2 is ghost cell, then cell_1 cannot be ghost
            disn(iface) = ABS((ghost_xc(cell_2 - ncells) - xc(cell_1))*vecfx(iface) + &
                            (ghost_yc(cell_2 - ncells) - yc(cell_1))*vecfy(iface) + &
                            (ghost_zc(cell_2 - ncells) - zc(cell_1))*vecfz(iface))
          ELSE
            disn(iface) = ABS((xc(cell_2) - xc(cell_1))*vecfx(iface) + &
                            (yc(cell_2) - yc(cell_1))*vecfy(iface) + &
                            (zc(cell_2) - zc(cell_1))*vecfz(iface))
          ENDIF
        END IF
    END DO
  ELSE
    DO iface=1,nfaces
        cell_1 =  lfc(iface,1)
        cell_2 =  lfc(iface,2)  
        IF(cell_1 == cell_2) THEN
          !...boundary face cannot be associated with a ghost cell
          disn(iface) = ABS((xf(iface) - xc(cell_1))*vecfx(iface) + &
                            (yf(iface) - yc(cell_1))*vecfy(iface))
        ELSE
          IF (cell_1 > ncells) THEN
            !...cell_1 is ghost cell, then cell_2 cannot be ghost
            disn(iface) = ABS((xc(cell_2) - ghost_xc(cell_1 - ncells))*vecfx(iface) + &
                            (yc(cell_2) - ghost_xc(cell_1 - ncells))*vecfy(iface))
          ELSEIF (cell_2 > ncells) THEN
            !...cell_2 is ghost cell, then cell_1 cannot be ghost
            disn(iface) = ABS((ghost_xc(cell_2 - ncells) - xc(cell_1))*vecfx(iface) + &
                            (ghost_yc(cell_2 - ncells) - yc(cell_1))*vecfy(iface))
          ELSE
            disn(iface) = ABS((xc(cell_2) - xc(cell_1))*vecfx(iface) + &
                            (yc(cell_2) - yc(cell_1))*vecfy(iface))
          ENDIF
        END IF
    END DO
  ENDIF
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished calculating disn"
  
  ! Finding orientation of normal at face relative to CV
  !...i.e. normdir(i,j) = 1 if face j is defined with cell_1 is i, otherwise normdir(i,j) = -1
  ALLOCATE (normdir(ncells,nf_max), stat=ierr)
  normdir(:,:) = 1
  DO i = 1,ncells
    DO j = 1,nfcell(i)
      iface = lcf(i,j)
      IF (lfc(iface,1) == i) THEN
        normdir(i,j) = 1   ! Outward already
      ELSE
        normdir(i,j) = -1  ! Flip to make outward
      END IF
    END DO
  END DO
  
  !CHANGED THE NORMDIR temporarily

  ! Compute Interpolant from cell to face
  ALLOCATE(wcf(nfaces), stat=ierr)
  DO iface = 1,nfaces
      cell_1 = lfc(iface,1)
      cell_2 = lfc(iface,2)
      IF (cell_1 > ncells) THEN
        !...cell_1 is ghost, then cell_2 cannot be ghost
        xtemp1 = SQRT((xf(iface) - ghost_xc(cell_1 - ncells))**2.0 + (yf(iface) - ghost_yc(cell_1 - ncells))**2.0)
        IF(threed) xtemp1 = SQRT((xf(iface) - ghost_xc(cell_1 - ncells))**2.0 + &
                              (yf(iface) - ghost_yc(cell_1 - ncells))**2.0 + &
                              (zf(iface) - ghost_zc(cell_1 - ncells))**2.0)
      ELSE
        xtemp1 = SQRT((xf(iface) - xc(cell_1))**2.0 + (yf(iface) - yc(cell_1))**2.0)
        IF(threed) xtemp1 = SQRT((xf(iface) - xc(cell_1))**2.0 + &
                                (yf(iface) - yc(cell_1))**2.0 + &
                                (zf(iface) - zc(cell_1))**2.0)
      ENDIF
      IF (cell_2 > ncells) THEN
        !...cell_2 is ghost, then cell_1 cannot be ghost
        xtemp2 = SQRT((xf(iface) - ghost_xc(cell_2 - ncells))**2.0 + (yf(iface) - ghost_yc(cell_2 - ncells))**2.0)
        IF(threed) xtemp2 = SQRT((xf(iface) - ghost_xc(cell_2 - ncells))**2.0 + &
                              (yf(iface) - ghost_yc(cell_2 - ncells))**2.0 + &
                              (zf(iface) - ghost_zc(cell_2 - ncells))**2.0)
      ELSE
        xtemp2 = SQRT((xf(iface) - xc(cell_2))**2.0 + (yf(iface) - yc(cell_2))**2.0)
        IF(threed) xtemp2 = SQRT((xf(iface) - xc(cell_2))**2.0 + &
                                (yf(iface) - yc(cell_2))**2.0 + &
                                (zf(iface) - zc(cell_2))**2.0)
      ENDIF
      wcf(iface) = xtemp2/(xtemp1 + xtemp2)
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished calculating wcf"
  
  ! Compute interpolant from cell to vertex
  !...as vertex weights are only for post-processing, I leave wcv & wbcv as is, will consider later...
  ALLOCATE(wcv(ncells,nv_max1), stat=ierr) 
  ALLOCATE(sumweight(nnodes), stat=ierr)
  ALLOCATE(bnode(nnodes), stat=ierr)

  sumweight(:) = zero
  wcv(:,:) = zero
  DO i=1,ncells
      DO j=1,nvcell(i)
        inode = lcv(i,j)
        disntemp = SQRT((xv(inode) - xc(i))**2.0 + (yv(inode) - yc(i))**2.0)
        IF(threed) disntemp = SQRT((xv(inode) - xc(i))**2.0 + &
                                    (yv(inode) - yc(i))**2.0 + &
                                    (zv(inode) - zc(i))**2.0)
        sumweight(inode) = sumweight(inode) + one/disntemp
        wcv(i,j) = one/disntemp
      END DO
  END DO

  DO i=1,ncells
      DO j=1,nvcell(i)
        inode = lcv(i,j)
        wcv(i,j) = wcv(i,j)/sumweight(inode)
      END DO
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished calculating wcv"
  DEALLOCATE(sumweight)

  ! Compute interpolant from boundary face to vertex
  ALLOCATE(wbfv(nbcfaces,nv_max2), stat=ierr) 
  ALLOCATE(sumwtf(nnodes), stat=ierr)
  
  sumwtf(:) = zero
  wbfv(:,:) = zero
  DO i=1,nbcfaces
    !IF (bctype(i) == CONJ) CYCLE 
    ifc = bf_to_f(i)
    DO j=1,nvface(ifc)
      inode = lfv(ifc,j)
      disntemp = SQRT((xv(inode) - xf(ifc))**2.0 + (yv(inode) - yf(ifc))**2.0)
      IF(threed) disntemp = SQRT((xv(inode) - xf(ifc))**2.0 + &
                                  (yv(inode) - yf(ifc))**2.0 + &
                                  (zv(inode) - zf(ifc))**2.0)
      sumwtf(inode) = sumwtf(inode) + one/disntemp
      wbfv(i,j) = one/disntemp
    END DO
  END DO

  DO i=1,nbcfaces
    !IF (bctype(i) == CONJ) CYCLE             !!! excluding conj face 30 july
    ifc = bf_to_f(i)
    DO j=1,nvface(ifc)
      inode = lfv(ifc,j)
      wbfv(i,j) = wbfv(i,j)/sumwtf(inode)
    END DO
  END DO
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') &
                      "Finished calculating weights"
  DEALLOCATE(sumwtf)

  ! Tagging the boundary nodes
  bnode(:) = 0
  DO iface = 1,nfaces
    ibc = f_to_bf(iface)
    ! IF (ibc .eq. 0) THEN
    !   WRITE(*,'(a,i0)') "iface= ",iface
    !   STOP
    ! ENDIF
    !IF ((bface(iface) == 1) .and. (bctype (ibc) /= CONJ)) THEN
    IF (bface(iface) == 1) THEN
      IF (bctype(ibc) /= CONJ) THEN
        DO j = 1,nvface(iface)
          inode = lfv(iface,j)
          bnode(inode) = 1
        END DO
      ENDIF
    END IF
  END DO
  
  ! Load temperature into boundary condition array
  ALLOCATE(temp_bc(nbcfaces))
  temp_bc(:) = -123456.d0
  DO ibc = 1,nbcfaces
    temp_bc(ibc) = phibc(ibc)
  ENDDO
  DEALLOCATE(phibc)

  !...should only rank 0 close the file?
  CLOSE(unit = io_ace)
  
  !...free memory of variables holding global mesh
  DEALLOCATE(nfcell_GB, nvcell_GB, owner_GB, gc_to_pgc, pgc_to_gc)
  DEALLOCATE(xc_GB, yc_GB, volcell_GB)
  DEALLOCATE(xf_GB, yf_GB)
  DEALLOCATE(vecfx_GB, vecfy_GB, areaf_GB)
  DEALLOCATE(xv_GB, yv_GB)

  DEALLOCATE(lcf_GB, lcv_GB, lfc_GB, lfv_GB, lcc_GB)
  DEALLOCATE(nvface_GB)
  DEALLOCATE(no_f_patch_GB)
  DEALLOCATE(phibc_GB)
  DEALLOCATE(bface_GB, f_to_bf_GB, bf_to_f_GB, bf_to_c_GB, bctype_GB)
  DEALLOCATE(cell_type_GB, material_type_GB)
  IF (threed) THEN
    DEALLOCATE(zc_GB)
    DEALLOCATE(zf_GB, vecfz_GB)
    DEALLOCATE(zv_GB)
  ENDIF

  DEALLOCATE(lf_to_gf, lv_to_gv, lbf_to_gbf)
  DEALLOCATE(gc_to_lc, gf_to_lf, gv_to_lv, gbf_to_lbf)
  DEALLOCATE(oc_to_gc)
  DEALLOCATE(ghost_cells_gc, ghost_cells_pgc)
  DEALLOCATE(ghost_xc, ghost_yc)
  IF (threed) THEN
    DEALLOCATE(ghost_zc)
  ENDIF
  
  IF ((debug_level > 0).and.(rank_cart == 0)) WRITE(io_debug,'(a)') "Exiting grid_and_boundary"
    
END SUBROUTINE grid_and_boundary