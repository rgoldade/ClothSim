add_executable(FEMHangingMesh FEMHangingMesh.cpp)

target_link_libraries(FEMHangingMesh
						PRIVATE
						ClothSolvers
						RenderTools
						Utilities)

file( RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )						

install(TARGETS FEMHangingMesh RUNTIME DESTINATION ${REL})
