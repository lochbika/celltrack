
set(target_name celltrack)
FILE(GLOB sourcefiles "*.f90")
add_executable(${target_name} ${sourcefiles})
install(TARGETS celltrack RUNTIME DESTINATION bin)

if(CDI_INCLUDE)
  include_directories(${CDI_INCLUDE})
else(CDI_INCLUDE)
  MESSAGE(STATUS "WARNING: No cdi includes are found.")
endif(CDI_INCLUDE)

if(CDI_LIB)
  target_link_libraries(celltrack ${CDI_LIB})
else(CDI_LIB)
  MESSAGE(STATUS "WARNING: No cdi lib was found.")
endif(CDI_LIB)
