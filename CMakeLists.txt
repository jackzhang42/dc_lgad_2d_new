# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(dc_lgad)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(dc_lgad dc_lgad.C)
target_link_libraries(dc_lgad Garfield::Garfield)

add_executable(pin pin.C)
target_link_libraries(pin Garfield::Garfield)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(dc_lgad OpenMP::OpenMP_CXX)
endif()

# ---Copy all files locally to the build directory-------------------------------
foreach(_file n3397_des.dat n3397_des.grd n3469_des.dat n3469_des.grd)
	configure_file(${_file} ${_file} COPYONLY)
endforeach()
