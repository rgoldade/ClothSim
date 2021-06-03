add_library(Utilities SimpleGeometryBuilder.cpp TriMesh.cpp)

target_link_libraries(Utilities PUBLIC RenderTools)

target_include_directories(Utilities PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
							  $<INSTALL_INTERFACE:include>)