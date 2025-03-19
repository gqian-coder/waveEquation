cmake_minimum_required(VERSION 3.19)

project(mgard_standard)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED) #PATHS /home/qian/Software/ADIOS2/install-adios/) 
find_package(mgard REQUIRED)


option(IOCOMP_ENABLE_HIP "Enable HIP compressors." ON)
set(CMAKE_CXX_FLAGS  "-fno-omit-frame-pointer -Wfatal-errors -Wall -Wextra")

#set (CMAKE_CXX_STANDARD 17)
#set (CMAKE_CXX_STANDARD_REQUIRED ON)

if (${IOCOMP_ENABLE_HIP})
    add_compile_definitions(ENABLE_HIP_COMPRESSOR=1)
endif()


add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -Wfatal-errors -Wall -Wextra")

add_executable(waveProp_twoStep waveProp_twoStep.cpp waveEquation.tpp waveInit.tpp obstacle.tpp)
target_link_libraries(waveProp_twoStep adios2::adios2 MPI::MPI_C mgard::mgard)

add_executable(waveProp_twoStep_3d waveProp_twoStep_3d.cpp waveEquation.tpp waveInit.tpp obstacle.tpp)
target_link_libraries(waveProp_twoStep_3d adios2::adios2 MPI::MPI_C mgard::mgard)

add_executable(waveEnergy waveEnergy.cpp)
target_link_libraries(waveEnergy adios2::adios2 MPI::MPI_C)

add_executable(waveEnergy_3d waveEnergy_3d.cpp)
target_link_libraries(waveEnergy_3d adios2::adios2 MPI::MPI_C)

add_executable(wavePotential wavePotential.cpp)
target_link_libraries(wavePotential adios2::adios2 MPI::MPI_C mgard::mgard)

