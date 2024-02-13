# Build the neptune executable based on the given source files
add_executable(neptune-valsent1b 	
							rdpoe.f90
							rdinp.f90
							valsent.f90)
add_executable(neptune-valsent1b-set 	
							rdpoe.f90
							rdinp.f90
							valsent_set.f90)

# Dynamic linking against our libraries
target_link_libraries(neptune-valsent1b 
						${LIBSLAM}
						neptune)
target_link_libraries(neptune-valsent1b-set 
						${LIBSLAM}
						neptune)

# Install the libraries into the lib directory
install(TARGETS neptune-valsent1b
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(TARGETS neptune-valsent1b-set
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# Install the Fortran module files into the include directory
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/NEPTUNE-valsent1b)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/NEPTUNE-valsent1b-set)