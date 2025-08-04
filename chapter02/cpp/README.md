testInstallationFC

A C++ project that automates dependency management for Boost and Eigen libraries while supporting modern C++20 standards. The project demonstrates efficient handling of dependencies by downloading and installing them if they are not already present.
Features

    Automatically downloads and installs Boost (1.87.0) to ~/boost_1_87_0 if not found.
    Automatically downloads and installs Eigen (3.4.0) to ~/eigen-3.4.0 if not found.
    Supports modern C++20 standards.
    Builds optimized binaries using the -Ofast..  compiler flag.
    Installs the compiled binaries to ~/bin by default.

Prerequisites

Before building the project, ensure the following are installed:

    CMake: Version 3.25 or higher.
    G++: A compiler that supports C++20 (e.g., GCC 10 or newer).
    GSL lib: install with e.g. : sudo apt install libgsl-dev 
    OMP lib: install with e.g. : sudo apt install libomp-dev
    Tools:
        tar: Required for extracting downloaded archives.
        Stable internet connection for downloading dependencies.

Getting Started
Clone the Repository

git clone <repository_url>
cd <repository_directory>

Build and Install

    Generate Build Files:

cmake -B build -DCMAKE_BUILD_TYPE=Release

Build the Project:

cmake --build build

Install the Binary:

    cmake --install build

    By default, the binary will be installed to /home/r/fortranWS/bin.

Customizing the Build
Custom Library Paths

    To use a custom Boost installation path:

cmake -B build -DBOOST_ROOT=/custom/path/to/boost

To use a custom Eigen installation path:

    cmake -B build -DEIGEN_INCLUDE_DIR=/custom/path/to/eigen

Change Binary Installation Path

    Modify the installation path during configuration:

    cmake -B build -DCMAKE_INSTALL_PREFIX=/custom/install/path

Specify Compiler

    To use a specific compiler, set it during configuration:

    cmake -B build -DCMAKE_CXX_COMPILER=/path/to/your/compiler

How It Works
Boost Library Handling

    Checks if BOOST_ROOT is defined and valid.
    If not, sets the default Boost path to ~/boost_1_87_0.
    Downloads and extracts Boost from https://archives.boost.io/release/1.87.0/source/ if the directory does not exist.

Eigen Library Handling

    Checks if EIGEN_INCLUDE_DIR is defined and valid.
    If not, sets the default Eigen path to ~/eigen-3.4.0.
    Downloads and extracts Eigen from https://gitlab.com/libeigen/eigen if the directory does not exist.


Checks if GSL_INCLUDE_DIR is defined and valid.
    If not: message : install  

Checks if OMP_DIR is defined and valid.
    If not: message : install  

Troubleshooting
Common Issues

    Boost or Eigen Not Found:
        Ensure your network connection is stable.
        Check for write permissions in the installation directories (~/boost_1_87_0 or ~/eigen-3.4.0).

    Download Failure:
        Verify the URLs for Boost and Eigen in the CMakeLists.txt file.
        Check proxy settings or firewall restrictions.

    CMake Errors:
        Ensure your CMake version meets the minimum requirement (3.25).
        Specify custom paths for dependencies during configuration.

Error Recovery

    Delete the build directory and rerun the configuration:

    rm -rf build
    cmake -B build

Extending the Project

    Add More Libraries:
        Follow the pattern used for Boost and Eigen to add more dependencies.
    Change Compiler Flags:
        Modify add_compile_options(-Ofast) in the CMakeLists.txt file for specific optimization needs.

FAQs
Q: What happens if Boost or Eigen is already installed?

If the specified directories (BOOST_ROOT or EIGEN_INCLUDE_DIR) exist, the script skips downloading and installing the libraries.
Q: Can I change the Boost or Eigen version?

Yes, modify the BOOST_VERSION or EIGEN_VERSION variables in CMakeLists.txt and update the corresponding URLs.
Q: Why is -Ofast used as a compiler flag?

The -Ofast flag enables aggressive optimizations for performance but may not strictly adhere to the C++ standard.

License
This project is licensed under the GPL 2.0 License. See the LICENSE file for details.
Acknowledgments

Special thanks to:

    Boost for providing robust C++ libraries.
    Eigen for offering efficient linear algebra tools.
    The open-source community for tools and support.

