cmake_minimum_required(VERSION 3.19)

project(mgard_standard)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED) #PATHS /home/qian/Software/ADIOS2/install-adios/) 
find_package(mgard REQUIRED)

find_package(zstd REQUIRED /lustre/orion/proj-shared/csc143/jieyang/MGARD/install-hip-frontier/lib/cmake/zstd/) 
set(ZSTD_INCLUDES "/lustre/orion/proj-shared/csc143/jieyang/MGARD/install-hip-frontier/include/")
find_library(ZSTD_LIB zstd HINTS "/lustre/orion/proj-shared/csc143/jieyang/MGARD/install-hip-frontier/lib/")

find_package(SZ3 REQUIRED PATHS /lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/SZ3/install/usr/local/lib64/cmake/SZ3/)
set(SZ3_INCLUDES "/lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/SZ3/install/usr/local/include")
find_library(SZ3_LIB SZ3c REQUIRED HINTS "/lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/SZ3/install/usr/local/lib64/")

find_package(zfp REQUIRED PATHS /lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/zfp/install/usr/local/)
set(ZFP_INCLUDES "/lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/zfp/install/usr/local/include")
find_library(ZFP_LIB zfp HINTS "/lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/zfp/install/usr/local/lib64/")


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

# ZFP
add_executable(waveProp_twoStep_3d_zfp
               waveProp_twoStep_3d_zfp.cpp
               waveEquation.tpp
               waveInit.tpp 
               obstacle.tpp)
target_include_directories(waveProp_twoStep_3d_zfp PRIVATE ${ZFP_INCLUDES})
target_link_libraries(waveProp_twoStep_3d_zfp adios2::adios2 MPI::MPI_C ${ZFP_LIB})

add_executable(waveProp_twoStep_zfp
               waveProp_twoStep_zfp.cpp
               waveEquation.tpp
               waveInit.tpp
               obstacle.tpp)
target_include_directories(waveProp_twoStep_zfp PRIVATE ${ZFP_INCLUDES})
target_link_libraries(waveProp_twoStep_zfp adios2::adios2 MPI::MPI_C ${ZFP_LIB})

# SZ
add_executable(waveProp_twoStep_3d_sz
               waveProp_twoStep_3d_sz.cpp
               waveEquation.tpp
               waveInit.tpp
               obstacle.tpp)
target_include_directories(waveProp_twoStep_3d_sz PRIVATE ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(waveProp_twoStep_3d_sz PRIVATE adios2::adios2 MPI::MPI_C ${ZSTD_LIB} ${SZ3_LIB})

add_executable(waveProp_twoStep_sz
               waveProp_twoStep_sz.cpp
               waveEquation.tpp
               waveInit.tpp
               obstacle.tpp)
target_include_directories(waveProp_twoStep_sz PRIVATE ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(waveProp_twoStep_sz PRIVATE adios2::adios2 MPI::MPI_C ${ZSTD_LIB} ${SZ3_LIB})

