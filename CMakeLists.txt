cmake_minimum_required(VERSION 3.12)

project(nays2dh Fortran)

include(GNUInstallDirs)

# set CMAKE_PREFIX_PATH to the directory containing iriclibConfig.cmake
# ie -DCMAKE_PREFIX_PATH=d:/Downloads/iricsdk/iriclib-4.0.x/lib/cmake/iriclib
find_package(iriclib 4 REQUIRED)

if (NOT WIN32)
  # Linux requires hdf5 and Poco to be located to link
  # set CMAKE_PREFIX_PATH to the directories containing iriclibConfig.cmake, hdf5-config.cmake, and PocoFoundationConfig.cmake
  # -DCMAKE_PREFIX_PATH=/iricdev-2019/lib/install/hdf5-1.8.21/share/cmake/hdf5;/iricdev-2019/lib/install/iriclib-4.0.14/lib/cmake;/iricdev-2019/lib/install/poco-1.9.4/lib/cmake/Poco
  find_package(hdf5 1.8 REQUIRED)
  find_package(PocoFoundation REQUIRED)
endif()

# Change this to ON if solver requires OpenMP
option(IRIC_SOLVER_USE_OPENMP "Solver requires OpenMP" ON)
if(IRIC_SOLVER_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# No need to add iric.f90 it will be added automatically by target_link_libraries
add_executable(nays2dh
  src/Nays2DH.f90
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if (WIN32)
    target_compile_options(nays2dh PRIVATE /nostandard-realloc-lhs)
  else()
    target_compile_options(nays2dh PRIVATE -nostandard-realloc-lhs)
  endif()
endif()

target_link_libraries(nays2dh PRIVATE iriclib::iriclib)

if(IRIC_SOLVER_USE_OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(nays2dh PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(TARGETS nays2dh DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(FILES nays2dh64/definition.xml DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES nays2dh64/translation_ja_JP.ts DESTINATION "${CMAKE_INSTALL_BINDIR}")

# include(CTest)

# add_subdirectory(tests)
