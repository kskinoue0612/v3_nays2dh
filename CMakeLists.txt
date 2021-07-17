cmake_minimum_required(VERSION 3.12)

project(<BASENAME> Fortran)

include(GNUInstallDirs)

# set CMAKE_PREFIX_PATH to the directory containing iriclibConfig.cmake
# ie -DCMAKE_PREFIX_PATH=d:/Downloads/iricsdk/iriclib-4.0.x/lib/cmake/iriclib
find_package(iriclib 4 REQUIRED)

# Change this to ON if solver requires OpenMP
option(IRIC_SOLVER_USE_OPENMP "Solver requires OpenMP" OFF)
if(IRIC_SOLVER_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# No need to add iric.f90 it will be added automatically by target_link_libraries
add_executable(<BASENAME>
  <ADD SOURCE FILES HERE>
)

target_link_libraries(<BASENAME> PRIVATE iriclib::iriclib)

if(IRIC_SOLVER_USE_OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(<BASENAME> PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(TARGETS <BASENAME> DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(FILES <LOCATION OF definition.xml> DESTINATION "${CMAKE_INSTALL_BINDIR}")

# include(CTest)

# add_subdirectory(tests)
