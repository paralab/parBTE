cmake_minimum_required(VERSION 3.15.3)
project(BTE)

enable_language(Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)

set(CMAKE_CXX_STANDARD 14)
FIND_PACKAGE(MPI REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)

if(OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release " FORCE)
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

set(SOURCE_FILES
    src/variables_module.f90
    src/precisions_module.f90
    src/grid_module.f90
    src/constants_module.f90
    src/petsc_solver_mod.f90        #... parent module of all 3 methods: matrix-assembled, aMat, aMat batch
    src/petsc_matAssembled_mod.f90  #... module for matrix-assembled method
    src/amat_parent_mod.f90         #... parent module of 2 methods: aMat, aMat batch
    src/amat_noBatched_mod.f90      #... module of aMat method
    src/amat_batched_mod.f90        #... module of aMat batch method
    src/aluminiumtemp.f90
    src/avgvelocity.f90
    src/epcalc.f90
    src/gauss_quad.f90
    src/grid_and_boundary.f90
    src/initialize.f90
    src/intergratedintensity.f90
    src/iwallcal.f90
    src/polar.f90
    src/postprocess_vtk.f90
    src/ppnodal.f90
    src/reflection.f90
    src/relaxtime_broido.f90
    src/siliconthermalcond.f90
    src/snweights.f90
    src/solver.f90
    src/tempcalc_notbc.f90
    src/tri.f90
    src/update.f90
    src/variabledt.f90
    src/utilities.f90
    src/bsnord_swa_mod.f90          #... use petsc matrix-assembled method (USE petsc_matAssembled)
    src/bsnord_swa_mod_aMat.f90     #... use aMat (no batch) method (USE amat_noBatched)
    src/bsnord_swa_mod_aMat_bch.f90 #... use aMat with batching method (USE amat_batched)
    src/profiler_module.f90
    src/avg_ttop_calc.f90)

option(BUILD_WITH_PETSC "Build code with the petsc" ON)

option(USE_PETSC_MATASSEMBLED "Build code using PETSc matrix-assembled method" OFF)
option(USE_AMAT_NOBATCHED "Build code using aMat noBatched method" ON)
option(USE_AMAT_BATCHED "Build code using aMat batched method" OFF)

if(BUILD_WITH_PETSC)
    list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
    find_package(PETSc REQUIRED)
    add_definitions(-DBUILD_WITH_PETSC)
endif(BUILD_WITH_PETSC)

if(USE_PETSC_MATASSEMBLED)
    add_definitions(-DUSE_PETSC_MATASSEMBLED)
endif(USE_PETSC_MATASSEMBLED)

if(USE_AMAT_NOBATCHED)
    add_definitions(-DUSE_AMAT_NOBATCHED)
endif(USE_AMAT_NOBATCHED)

if(USE_AMAT_BATCHED)
    add_definitions(-DUSE_AMAT_BATCHED)
endif(USE_AMAT_BATCHED)

add_executable(tdtr_no_ksub_broi ${SOURCE_FILES})

target_include_directories(tdtr_no_ksub_broi PUBLIC include)
target_include_directories(tdtr_no_ksub_broi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(tdtr_no_ksub_broi ${MPI_Fortran_LIBRARIES} m)

if(BUILD_WITH_PETSC)
    target_include_directories(tdtr_no_ksub_broi PUBLIC ${PETSC_DIR}/include)
    target_include_directories(tdtr_no_ksub_broi PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(tdtr_no_ksub_broi ${PETSC_LIBRARIES})
endif(BUILD_WITH_PETSC)