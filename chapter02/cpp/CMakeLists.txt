cmake_minimum_required(VERSION 3.25)

# Define the project name and version
project(unitTestGSL VERSION 1.0)

# Set the C++ standard and explicitly specify the compiler
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER /usr/bin/g++)

# Set installation prefix to ~/bin
set(INSTALL_PREFIX "$ENV{HOME}/bin")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE PATH "Default install path" FORCE)

# Include directories for project-specific headers
set(INCLUDE_DIR "./include")

# Compiler optimization flags
add_compile_options(-Ofast -march=native -Wall -fopenmp)

# --- Check for Required Libraries ---

# Ensure OpenBLAS and LAPACKE are installed
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(FATAL_ERROR 
        "OpenBLAS or LAPACKE library not found. Please install with:\n"
        "sudo apt-get install libopenblas-dev liblapacke-dev"
    )
endif()

# Function to install missing libraries
function(install_library LIB_NAME LIB_VERSION LIB_URL LIB_INSTALL_DIR_VAR)
    set(DEFAULT_INSTALL_DIR "$ENV{HOME}/${LIB_NAME}_${LIB_VERSION}")
    set(${LIB_INSTALL_DIR_VAR} "${DEFAULT_INSTALL_DIR}" PARENT_SCOPE)

    if(NOT EXISTS "${DEFAULT_INSTALL_DIR}")
        message(STATUS "${LIB_NAME} not found in ${DEFAULT_INSTALL_DIR}, downloading...")

        file(MAKE_DIRECTORY "${DEFAULT_INSTALL_DIR}")
        set(DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/${LIB_NAME}_download")
        file(MAKE_DIRECTORY "${DOWNLOAD_DIR}")

        set(ARCHIVE_NAME "${LIB_NAME}-${LIB_VERSION}.tar.gz")
        set(ARCHIVE_PATH "${DOWNLOAD_DIR}/${ARCHIVE_NAME}")
        if(NOT EXISTS "${ARCHIVE_PATH}")
            file(DOWNLOAD "${LIB_URL}" "${ARCHIVE_PATH}" SHOW_PROGRESS STATUS dl_status)
            list(GET dl_status 0 dl_result)
            if(NOT "${dl_result}" STREQUAL "0")
                message(FATAL_ERROR "Failed to download ${LIB_NAME}. Check your network or URL.")
            endif()
        endif()

        execute_process(
            COMMAND tar -xzf "${ARCHIVE_PATH}" -C "${DEFAULT_INSTALL_DIR}" --strip-components=1
            RESULT_VARIABLE extract_result
        )
        if(NOT "${extract_result}" STREQUAL "0")
            message(FATAL_ERROR "Failed to extract ${LIB_NAME}.")
        endif()
    else()
        message(STATUS "${LIB_NAME} is already installed in ${DEFAULT_INSTALL_DIR}.")
    endif()
endfunction()

# --- Library Installations --- 
# Boost Configuration
install_library(
    "boost" "1.88.0"
    "https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz"
    BOOST_INCLUDE_DIR
)

# Eigen Configuration
install_library(
    "eigen" "3.4.0"
    "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
    EIGEN_INCLUDE_DIR
)

# Add library include paths
include_directories(${INCLUDE_DIR} ${BOOST_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})

# --- Find and Link Libraries ---

# GSL (Now Optional)
find_package(GSL)
if(GSL_FOUND)
    message(STATUS "GSL library found: ${GSL_INCLUDE_DIRS}")
    set(USE_GSL TRUE)
else()
    message(WARNING "GSL library not found. Some features will be disabled. Install with: sudo apt-get install gsl-bin libgsl0-dev")
    set(USE_GSL FALSE)
endif()

# OpenMP (Required)
find_package(OpenMP REQUIRED)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "OpenMP library not found. Ensure it's available in your compiler.")
endif()

# --- Source Files ---
set(SOURCES
   # ompGammaBoostD.cpp
   # ompGammaBoostL1024.cpp
   # ompGammaBoostL125.cpp
   # ompGammaBoostL1331.cpp
    ompGammaGslD.cpp
   # ompGammaLanD.cpp
    ompGammaLibD.cpp
    ompGammaLibF.cpp
    ompGammaLibL.cpp
    ompGammaSerD.cpp
    ompGammaSerF.cpp
    ompGammaSerL27.cpp
    ompGammaSerL30.cpp
    ompGammaSerL44.cpp
)

# Loop through the source files and create executables
foreach(SOURCE_FILE ${SOURCES})
    get_filename_component(TARGET_NAME ${SOURCE_FILE} NAME_WE)

    add_executable(${TARGET_NAME} ${SOURCE_FILE})

    # Link required libraries
    target_link_libraries(${TARGET_NAME} OpenMP::OpenMP_CXX BLAS::BLAS LAPACK::LAPACK)
    
    # Link GSL only if found
    if(USE_GSL)
        target_link_libraries(${TARGET_NAME} GSL::gsl GSL::gslcblas)
    endif()

    install(TARGETS ${TARGET_NAME} DESTINATION ${INSTALL_PREFIX})
endforeach()

message(STATUS "Configuration complete. Use 'make' to build the project.")
