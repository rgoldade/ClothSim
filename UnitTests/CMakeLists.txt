if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/External/googletest" googletest)

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} PRIVATE gmock gtest gtest_main ClothSolvers Utilities)
    set_target_properties(${TESTNAME} PROPERTIES FOLDER UnitTests)
endmacro()

package_add_test(ClothSimTests TK20EnergyTests.cpp TriMeshTests.cpp)