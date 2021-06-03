add_library(ClothSolvers ExplicitFEMSolver.cpp FEMClothMesh.cpp TK20Energy.cpp)

target_link_libraries(ClothSolvers PUBLIC RenderTools Utilities)

target_include_directories(ClothSolvers PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
							  $<INSTALL_INTERFACE:include>)

target_include_directories(ClothSolvers SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/External/autodiff")