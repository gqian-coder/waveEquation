cmake_minimum_required(VERSION 3.19)

project(mgard_standard)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED PATHS /home/qian/Software/ADIOS2/install-adios/) 
find_package(mgard REQUIRED)

find_package(CUDAToolkit REQUIRED)
set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -w")

add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -Wfatal-errors -Wall -Wextra")

#add_executable(waveProp waveProp.cpp waveEquation.tpp)
#target_link_libraries(waveProp adios2::adios2 MPI::MPI_C mgard::mgard)

add_executable(waveProp_twoStep waveProp_twoStep.cpp waveEquation.tpp waveInit.tpp obstacle.tpp)
target_link_libraries(waveProp_twoStep adios2::adios2 MPI::MPI_C mgard::mgard)

add_executable(waveProp_twoStep_3d waveProp_twoStep_3d.cpp waveEquation.tpp waveInit.tpp obstacle.tpp)
target_link_libraries(waveProp_twoStep_3d adios2::adios2 MPI::MPI_C mgard::mgard)

add_executable(waveEnergy waveEnergy.cpp)
target_link_libraries(waveEnergy adios2::adios2 MPI::MPI_C)

#add_executable(waveProp_twoStep_m1 waveProp_twoStep_m1.cpp waveEquation.tpp waveInit.tpp obstacle.tpp)
#target_link_libraries(waveProp_twoStep_m1 adios2::adios2 MPI::MPI_C mgard::mgard)

#add_executable(waveProp_mgr waveProp_mgr.cpp waveEquation.tpp)
#target_link_libraries(waveProp_mgr adios2::adios2 MPI::MPI_C mgard::mgard)

