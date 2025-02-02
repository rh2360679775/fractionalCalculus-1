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

# Function to handle library installation
function(install_library LIB_NAME LIB_VERSION LIB_URL LIB_INSTALL_DIR_VAR)
    set(DEFAULT_INSTALL_DIR "$ENV{HOME}/${LIB_NAME}_${LIB_VERSION}")
    set(${LIB_INSTALL_DIR_VAR} "${DEFAULT_INSTALL_DIR}" PARENT_SCOPE)

    if(NOT EXISTS "${DEFAULT_INSTALL_DIR}")
        message(STATUS "${LIB_NAME} not found in ${DEFAULT_INSTALL_DIR}, downloading...")

        # Create installation and download directories
        file(MAKE_DIRECTORY "${DEFAULT_INSTALL_DIR}")
        set(DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/${LIB_NAME}_download")
        file(MAKE_DIRECTORY "${DOWNLOAD_DIR}")

        # Download library if not already downloaded
        set(ARCHIVE_NAME "${LIB_NAME}-${LIB_VERSION}.tar.gz")
        set(ARCHIVE_PATH "${DOWNLOAD_DIR}/${ARCHIVE_NAME}")
        if(NOT EXISTS "${ARCHIVE_PATH}")
            file(DOWNLOAD "${LIB_URL}" "${ARCHIVE_PATH}" SHOW_PROGRESS STATUS dl_status)
            list(GET dl_status 0 dl_result)
            if(NOT "${dl_result}" STREQUAL "0")
                message(FATAL_ERROR "Failed to download ${LIB_NAME}. Check your network or URL.")
            endif()
        endif()

        # Extract archive
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
    "boost" "1.86.0"
    "https://boostorg.jfrog.io/artifactory/main/release/1.86.0/source/boost_1_86_0.tar.gz"
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
find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL library not found. Install with: sudo apt-get install gsl-bin libgsl0-dev")
endif()

find_package(OpenMP REQUIRED)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "OpenMP library not found. Ensure it's available in your compiler.")
endif()

# --- Source Files ---
# List of source files
set(SOURCES
#ompGammaBoostD.cpp
#ompGammaBoostL.cpp
#ompGammaBoostL1024.cpp
#ompGammaBoostL125.cpp
#ompGammaBoostL1331.cpp
#ompGammaBoostL64.cpp
ompGammaGslD.cpp
ompGammaLanD.cpp
ompGammaLibD.cpp
ompGammaLibF.cpp
ompGammaLibL.cpp
ompGammaSerD.cpp
ompGammaSerF.cpp
ompGammaSerL27.cpp
ompGammaSerL30.cpp
ompGammaSerL44.cpp
testGlobalPadeMLF.cpp
)

# Loop through the source files and create executables
foreach(SOURCE_FILE ${SOURCES})
    # Extract the filename without extension to use as target name
    get_filename_component(TARGET_NAME ${SOURCE_FILE} NAME_WE)

    # Create an executable
    add_executable(${TARGET_NAME} ${SOURCE_FILE})

    # Link required libraries
    target_link_libraries(${TARGET_NAME} GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX )

    # Define installation rule for the executable
    install(TARGETS ${TARGET_NAME} DESTINATION ${INSTALL_PREFIX})
endforeach()

# --- GPL-2.0 License Information ---
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
