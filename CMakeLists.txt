cmake_minimum_required (VERSION 2.8.11)
project (emp-ag2pc)
set(NAME "emp-ag2pc")
set(CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "" FORCE)

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)

include(${CMAKE_FOLDER}/cmake/common.cmake)
include(${CMAKE_FOLDER}/cmake/source_of_randomness.cmake)
include(${CMAKE_FOLDER}/cmake/threading.cmake)

FIND_PACKAGE(emp-ot REQUIRED)
INCLUDE_DIRECTORIES(${EMP-OT_INCLUDE_DIRS})

install(DIRECTORY emp-ag2pc DESTINATION include)
install(DIRECTORY cmake/ DESTINATION cmake)

# Test cases
macro (add_test _name)
	add_test_with_lib(${_name} ${EMP-OT_LIBRARIES})
endmacro()


#add_test(aes)
#add_test(sha1)
#add_test(sha256)
#add_test (simple_circuit)
#add_test(demo demo)
#add_test(demo_parallel demo_parallel)
add_test(demo_1worker_sort demo_1worker_sort)
add_test(demo_1worker_test demo_1worker_test)
add_test(demo_1worker_custom demo_1worker_custom)
add_test(demo_orig_test demo_orig_test)
add_test(demo_orig_custom demo_orig_custom)
#add_test(mpi_test mpi_test)


SET (Boost_ADDITIONAL_VERSIONS "1.48" "1.48.0" "1.49" "1.49.0" "1.50" "1.50.0" "1.56" "1.56.0" "1.60.0")
# FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
FIND_PACKAGE (Boost COMPONENTS "mpi" "serialization" REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#target_link_libraries(demo ${Boost_LIBRARIES})
#target_link_libraries(demo ${MPI_CXX_LIBRARIES})
#target_link_libraries(demo_parallel ${Boost_LIBRARIES})
#target_link_libraries(demo_parallel ${MPI_CXX_LIBRARIES})
target_link_libraries(demo_1worker_sort ${Boost_LIBRARIES})
target_link_libraries(demo_1worker_sort ${MPI_CXX_LIBRARIES})
target_link_libraries(demo_1worker_test ${Boost_LIBRARIES})
target_link_libraries(demo_1worker_test ${MPI_CXX_LIBRARIES})
target_link_libraries(demo_1worker_custom ${Boost_LIBRARIES})
target_link_libraries(demo_1worker_custom ${MPI_CXX_LIBRARIES})
target_link_libraries(demo_orig_test ${Boost_LIBRARIES})
target_link_libraries(demo_orig_test ${MPI_CXX_LIBRARIES})
target_link_libraries(demo_orig_custom ${Boost_LIBRARIES})
target_link_libraries(demo_orig_custom ${MPI_CXX_LIBRARIES})
#target_link_libraries(mpi_test ${Boost_LIBRARIES})
#target_link_libraries(mpi_test ${MPI_CXX_LIBRARIES})
