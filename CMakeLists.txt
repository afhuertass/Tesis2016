CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(foo)
FIND_PACKAGE(ArrayFire REQUIRED)
FIND_PACKAGE(CUDA QUIET)
INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS}  )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")




IF (${CUDA_FOUND})
    IF(${ArrayFire_CUDA_FOUND})  # variable defined by FIND(ArrayFire ...)
        FIND_LIBRARY( CUDA_NVVM_LIBRARY
          NAMES "nvvm"
          PATH_SUFFIXES "nvvm/lib64" "nvvm/lib"
          PATHS ${CUDA_TOOLKIT_ROOT_DIR}
          DOC "CUDA NVVM Library"
          )
        MESSAGE(STATUS "EXAMPLES: CUDA backend is ON.")
	
	set( SDIR src )
	set( sources  ${SDIR}/main.cpp ${SDIR}/ForceModule.cpp  ${SDIR}/LBD3Q15.cpp ${SDIR}/ObstacleLoader.cpp ${SDIR}/Vector.cpp ${SDIR}/ELBD3Q15.cpp  )
	
	set( IDIR inc )
	set( headers ${IDIR}/ForceModule.h ${IDIR}/LBD3Q15.h ${IDIR}/ObstacleLoader.h ${IDIR}/Vector.h ${IDIR}/ELB3DQ15.h )
	
	
	include_directories( ${IDIR}/ )
	
	get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
	foreach(dir ${dirs})
	  message(STATUS "dir='${dir}'")
	endforeach()
	
	ADD_EXECUTABLE(lattice ${sources} ${headers} )
	TARGET_LINK_LIBRARIES(lattice ${ArrayFire_CUDA_LIBRARIES} "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES};${CUDA_cusolver_LIBRARY};${CUDA_CUFFT_LIBRARIES};${CUDA_NVVM_LIBRARY};${CUDA_CUDA_LIBRARY}" )
	
	ENDIF()
	ENDIF()


#ADD_EXECUTABLE(lattice 2LBD3Q15.cpp 2latticed3q15.h )
#TARGET_LINK_LIBRARIES(lattice  ${ArrayFire_CUDA_LIBRARIES} ${EXTRA_LIBS} )
