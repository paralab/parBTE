#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscmat.h>
SUBROUTINE Initialize
  USE petscmat
  USE precisions 
  USE variables
  USE grid, only : ncells, nnodes, nbcfaces, nfaces, temp_bc
  USE constants, only : zero, pi, four, one
  
  !USE aMat
!#ifdef USE_PETSC_MATASSEMBLED
  USE petsc_matAssembled
!#elif USE_AMAT_NOBATCHED
  USE amat_noBatched
!#elif USE_AMAT_BATCHED
  USE amat_batched
!#endif
  
  Implicit none
  
  integer(int_p) :: dim
  
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting initialize" 

  !nbands = nla + nta + nlo + nto   
  dim = ncells

  !IF (matMethod == 0) THEN
#if defined USE_PETSC_MATASSEMBLED
  CALL allocate_petsc_rhs_n_sol !...so far, aMat holds these vectors, later will move outside of aMat module
  CALL allocate_petsc_mat
  !CALL print_vectype(mp_rhs)
  !CALL print_vectype(mp_sol)
  !CALL print_mattype(mp_mat)
  !ELSEIF (matMethod == 1) THEN
#elif USE_AMAT_NOBATCHED
  CALL init_aMat                !...allocate variables used in aMat approach, and identify dep/indep cells
  CALL allocate_petsc_rhs_n_sol !...so far, aMat holds these vectors, later will move outside of aMat module
  CALL allocate_petsc_matshell  !...create matrix shell
  CALL buildScatterMap          !...build scatter map
  !ELSEIF (matMethod == 2) THEN
#elif USE_AMAT_BATCHED
  CALL init_aMat_bch(batchSize)
  CALL allocate_amat_rhs_n_sol_bch
  CALL allocate_petsc_matshell_bch
  CALL buildScatterMap_bch
  ! ELSEIF (matMethod == 3) THEN
  !   CALL allocate_petsc_rhs_n_sol_bch(batchSize)
  !   CALL allocate_petsc_mat_bch
  !   CALL build_maps
  !ELSE
#else
  WRITE(*,'(a)') "Initialize: wrong matrix method"
  STOP
  !ENDIF
#endif

  !...create ksp associated with assembled matrix or matrix shell
#if defined USE_PETSC_MATASSEMBLED
  CALL create_petsc_ksp
#elif USE_AMAT_NOBATCHED
  CALL create_petsc_ksp
#elif USE_AMAT_BATCHED
  CALL create_petsc_ksp_bch
#endif

  !... create petsc matrix, number of nonzeros per row is 7 (1 diagonal + 6 faces)
  ! CALL MatCreate(comm_cell, p_mat, p_ierr); CHKERRA(p_ierr)
  ! CALL MatSetsizes(p_mat, ncells, ncells, PETSC_DECIDE, PETSC_DECIDE, p_ierr); CHKERRA(p_ierr)
  ! IF (numprocs_cell > 1) THEN
  !   CALL MatSetType(p_mat, MATMPIAIJ, p_ierr); CHKERRA(p_ierr)
  !   CALL MatMPIAIJSetPreallocation(p_mat, 5, PETSC_NULL_INTEGER, 5, PETSC_NULL_INTEGER, p_ierr); CHKERRA(p_ierr) !...for tet element, max 4 faces per cell, plus itself
  ! ELSE
  !   CALL MatSetType(p_mat, MATSEQAIJ, p_ierr); CHKERRA(p_ierr)
  !   CALL MatSeqAIJSetPreallocation(p_mat, 5, PETSC_NULL_INTEGER, p_ierr); CHKERRA(p_ierr)
  ! ENDIF
  
  
  !... create petsc vector for RHS and solution
  ! CALL VecCreate(comm_cell, p_rhs, p_ierr); CHKERRA(p_ierr)
  ! IF (numprocs_cell > 1) THEN
  !   CALL VecSetType(p_rhs, VECMPI, p_ierr); CHKERRA(p_ierr)
  ! ELSE
  !   CALL VecSetType(p_rhs, VECSEQ, p_ierr); CHKERRA(p_ierr)
  ! ENDIF
  ! CALL VecSetSizes(p_rhs, ncells, PETSC_DECIDE, p_ierr); CHKERRA(p_ierr)
  ! CALL VecDuplicate(p_rhs, p_sol, p_ierr); CHKERRA(p_ierr)
  

  !... create ksp context
  ! CALL KSPCreate(comm_cell, p_ksp, p_ierr); CHKERRA(p_ierr)

  !... set matrix associated with the linear system
  ! CALL KSPSetOperators(p_ksp, p_mat, p_mat, p_ierr)
  ! CALL KSPSetType(p_ksp, KSPGMRES, p_ierr)
  ! CALL KSPSetFromOptions(p_ksp, p_ierr)
  ! CALL KSPSetTolerances(p_ksp, tol_inner, PETSC_DEFAULT_REAL, &
  !      PETSC_DEFAULT_REAL, iter_gmres, p_ierr) !...rtol, atol, dtol, max_iterations
  ! CALL KSPGetPC(p_ksp, p_pc, p_ierr);
  ! CALL PCSetType(p_pc, PCLU, p_ierr)
  ! CALL PCSetFromOptions(p_pc, p_ierr)
  


  ALLOCATE(tnot(ncells), tnotv(nnodes), tprime(ncells), tbprime(nbcfaces), &
    tnotone(ncells), temp_bcone(nbcfaces), ttop(nbcfaces), tbot(nbcfaces), ttopone(nbcfaces))
  ALLOCATE(gc(ncells, np_max, nbands), gcone(ncells, np_max, nbands))
  ALLOCATE(s(ncells))
  !...ALLOCATE(resid(dim))
  !...AlLOCATE(r(dim), soln(dim))
  !...ALLOCATE(aps(dim), anbs(dim,6), scs(dim,4)) !...scs_tran is not used at all, also for tet cell we only need anbs(dim,4)
  !...ALLOCATE(am(dim*7), ja(dim*7), ia(dim+1))
  !...ALLOCATE(ju(dim*7), vv(dim,21), jlu(dim*7), jr(2*dim), alu(dim*7)) 
  ALLOCATE(intensity(ncells, ndir, (highband-lowband+1), np_max),&
    intensityone(ncells, ndir, (highband-lowband+1), np_max))      !! dont have np_max, assuming 2 types of polarization, SI,GE
  !ALLOCATE(polar(nbands))
  ALLOCATE(gpvel(nbands, np_max))            !! assuming la and ta phonons only, silicon and germaniun 
  ALLOCATE(centfreq(nbands))
  
  ALLOCATE(Omega(ndir))
  ALLOCATE(rmu(ndir), rxi(ndir), ret(ndir))
  ALLOCATE(inrmu(ndir), inrxi(ndir), inret(ndir))      
 
  ALLOCATE(flux(nbcfaces, np_max, nbands), flux_s(nbcfaces, np_max, nbands))
  ALLOCATE(jfacepos(nbcfaces, np_max, nbands), jfaceneg(nbcfaces, np_max, nbands))
  ALLOCATE(jfacepos_s(nbcfaces, np_max, nbands), jfaceneg_s(nbcfaces, np_max, nbands))		!!!!sendbuffer
  ALLOCATE(inot_c(ncells, nbands, np_max))
  ALLOCATE(no_polar(nbands))
  ALLOCATE(delta_band(nbands))
  ALLOCATE(band_low(nbands))
  !ALLOCATE(checkpolar(2,nbands,4)) 
  !ALLOCATE(mat_polar(2,4), wavenumber(nbands,4))
  ALLOCATE(wavenumber(nbands, np_max)) !...changed from 4 to np_max

  !ALLOCATE(sdotn_array(nfaces, ndir))
  !ALLOCATE(insdotn_array(nfaces, ndir))

  ALLOCATE(qgen(ncells)) 
  !!! interior flux
  !ALLOCATE(flux_crossfar(4,nbands), flux_tranfar(4,nbands))
  !ALLOCATE(flux_crossnear(4,nbands), flux_trannear(4,nbands))
  !ALLOCATE(flux_crossfar_s(4,nbands), flux_tranfar_s(4,nbands))
  !ALLOCATE(flux_crossnear_s(4,nbands), flux_trannear_s(4,nbands))

  ALLOCATE(k_sil(ncells))

  s(:) = zero   
  tnot(:) = initial_temp
  tnotone(:) = initial_temp
  tprime(:) = zero
  tbprime(:) = zero 
  gc(:,:,:) = zero
  gcone(:,:,:) = zero
  intensity(:,:,:,:) = zero
  intensityone(:,:,:,:) = zero
  flux(:,:,:) = zero
  jfacepos(:,:,:) = zero
  jfaceneg(:,:,:) = zero
  inot_c(:,:,:) = zero

  flux_s(:,:,:) = zero
  jfacepos_s(:,:,:) = zero
  jfaceneg_s(:,:,:) = zero
 
  ! flux_crossfar(:,:)= zero
  ! flux_tranfar(:,:)= zero
  ! flux_crossnear(:,:)= zero
  ! flux_trannear(:,:)= zero 
  ! flux_crossfar_s(:,:)= zero
  ! flux_tranfar_s(:,:)= zero
  ! flux_crossnear_s (:,:)= zero
  ! flux_trannear_s(:,:)= zero
 
  qgen(:) = zero

  k_sil(:) = zero
  ttop(:) = initial_temp
  ttopone(:) = initial_temp
  tbot(:) = zero
  
  gpvel(:,:) = zero     ! initializin group velocity to be zero.
  
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing intialize" 

END SUBROUTINE Initialize