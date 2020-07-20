# Build the neptune executable based on the given source files
add_executable(neptune-valsent1b 	
							rdpoe.f90
							rdinp.f90
							valsent.f90)

# Dynamic linking against our libraries
target_link_libraries(neptune-valsent1b 
						${LIBSLAM}
						neptune)

# Install the libraries into the lib directory
install(TARGETS neptune-valsent1b
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# Install the Fortran module files into the include directory
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/NEPTUNE-valsent1b)