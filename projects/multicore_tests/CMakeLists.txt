#=============================================================================
# CMake configuration file for projects requiring the Chrono MULTICORE module
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message("\n---- MULTICORE projects")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    test_MCORE_foam
    test_MCORE_objectdrop
    test_MCORE_soilbin
    test_MCORE_suspension
    test_MCORE_wheel
    test_MCORE_radImSchmutz
    test_MCORE_settling
)

set(DEMOS_OPENGL
    test_MCORE_active_box
)
#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode.

find_package(Chrono QUIET
             COMPONENTS Multicore
             OPTIONAL_COMPONENTS OpenGL
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for MULTICORE projects")
  return()
endif()

if (CHRONO_OPENGL_FOUND)
    set(DEMOS ${DEMOS} ${DEMOS_OPENGL})
endif()

#--------------------------------------------------------------
# Include paths and libraries
#--------------------------------------------------------------

# (A) Path to the Chrono include headers
# - If using an installed version of Chrono, this will be the path 
#   to the installed headers (the configuration headers are also
#   available there)
# - If using a build version of Chrono, this will contain both the
#   path to the Chrono sources and the path to the chrono BUILD tree
#   (the latter for the configuration headers)
# 
# (B) Path to the top of the source tree for this project
# - for access to utility headers

include_directories(
    ${CHRONO_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
)

#--------------------------------------------------------------
# Append to the parent's list of DLLs (and make it visible up)
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Multicore tests...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
  target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach()
