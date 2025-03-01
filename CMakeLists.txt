cmake_minimum_required(VERSION 3.6)
project(GEDc++)


set(CMAKE_CXX_STANDARD 17)
option(CXX "enable C++ compilation" ON)
if (CXX)
    enable_language(CXX)
endif ()


# a CMake module named "FindGUROBI.cmake" is available in cmake/modules/
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})


add_subdirectory(src/utils)
add_subdirectory(src/executables)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


if (DEFINED ENV{GEDLIB_ROOT})
    message(STATUS "GEDLIB_ROOT in ENV I see is: $ENV{GEDLIB_ROOT}")
endif ()
if (NOT DEFINED ENV{GEDLIB_ROOT})
    message(STATUS "!!! GEDLIB_ROOT not set, please add GEDLIB_ROOT=~/gedlib to your cmake settings (e.g. in CLion), or via -DGEDLIB_ROOT=~/gedlib !!!")
endif ()

if (DEFINED ENV{LIBLSAP_ROOT})
    message(STATUS "LIBLSAP_ROOT in ENV I see is: $ENV{LIBLSAP_ROOT}")
endif ()
if (NOT DEFINED ENV{LIBLSAP_ROOT})
    message(STATUS "!!! LIBLSAP_ROOT not set, please add LIBLSAP_ROOT=~/liblsap to your cmake settings (e.g. in CLion), or via -DLIBLSAP_ROOT=~/liblsap !!!")
endif ()


add_executable(gurobi_muta
        src/executables/gurobi_muta.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(gurobi_muta
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(gurobi_muta
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(gurobi_muta
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        ${SCIP_LIBRARIES}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)

add_executable(gurobi_aids
        src/executables/gurobi_aids.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(gurobi_aids
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(gurobi_aids
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(gurobi_aids
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        ${SCIP_LIBRARIES}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)


add_executable(gurobi_protein
        src/executables/gurobi_protein.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(gurobi_protein
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(gurobi_protein
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(gurobi_protein
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        ${SCIP_LIBRARIES}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)


add_executable(gurobi_imdb
        src/executables/gurobi_imdb.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(gurobi_imdb
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(gurobi_imdb
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(gurobi_imdb
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        ${SCIP_LIBRARIES}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)


add_executable(subgraph_pattern
        src/executables/subgraph_pattern.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(subgraph_pattern
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)

target_link_directories(subgraph_pattern
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)

target_link_libraries(subgraph_pattern
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        ${SCIP_LIBRARIES}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)


set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
install(TARGETS subgraph_pattern
        RUNTIME
        DESTINATION bin)
install(TARGETS gurobi_muta
        RUNTIME
        DESTINATION bin)
install(TARGETS gurobi_aids
        RUNTIME
        DESTINATION bin)
install(TARGETS gurobi_protein
        RUNTIME
        DESTINATION bin)
install(TARGETS gurobi_imdb
        RUNTIME
        DESTINATION bin)
# overwrite install() command with a dummy macro that is a nop
macro(install)
endmacro()

