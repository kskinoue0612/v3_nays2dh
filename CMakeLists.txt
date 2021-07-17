cmake_minimum_required(VERSION 3.12)

project(nays2dh Fortran)

include(GNUInstallDirs)

# set CMAKE_PREFIX_PATH to the directory containing iriclibConfig.cmake
# ie -DCMAKE_PREFIX_PATH=d:/Downloads/iricsdk/iriclib-4.0.x/lib/cmake/iriclib
find_package(iriclib 4 REQUIRED)

# Change this to ON if solver requires OpenMP
option(IRIC_SOLVER_USE_OPENMP "Solver requires OpenMP" ON)
if(IRIC_SOLVER_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# No need to add iric.f90 it will be added automatically by target_link_libraries
add_executable(nays2dh
  src/Nays2DH.f90
)

target_link_libraries(nays2dh PRIVATE iriclib::iriclib)

if(IRIC_SOLVER_USE_OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(nays2dh PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(TARGETS nays2dh DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(FILES nays2dh64/definition.xml DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES nays2dh64/translation_ja_JP.ts DESTINATION "${CMAKE_INSTALL_BINDIR}")

# include(CTest)

# add_subdirectory(tests)
