cmake_minimum_required(VERSION 3.31)

project(tope LANGUAGES C)

# Only debug and release:
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

if(UNIX)
  set(CMAKE_C_FLAGS "-fvisibility=hidden -fPIC -no-pie -Wall -Wextra -std=c99 -pedantic" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -D_DEBUG -g" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG -s" CACHE STRING "" FORCE)
endif()

if (MSVC)
  # Replace dynamic runtime library (/MD) with static runtime library (/MT):
  set(variables CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE)
  foreach(variable ${variables})
    string(REGEX REPLACE "/MD" "/MT" var "${${variable}}")
    set(${variable} ${var} CACHE STRING "" FORCE)
  endforeach()

  # Suppress some warnings:
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4996 /wd6255 /wd6011")

  # Missing function declaration is error:
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /we4013")

  # No manifest:
  set(CMAKE_EXE_LINKER_FLAGS "/MANIFEST:NO" CACHE STRING "" FORCE)  # no EXE manifest
  set(CMAKE_SHARED_LINKER_FLAGS "/MANIFEST:NO" CACHE STRING "" FORCE)  # no DLL manifest

  # Remove unused flags:
  unset(CMAKE_C_FLAGS_MINSIZEREL CACHE)
  unset(CMAKE_C_FLAGS_RELWITHDEBINFO CACHE)
  unset(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL CACHE)
  unset(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO CACHE)
  unset(CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL CACHE)
  unset(CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO CACHE)
  unset(CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL CACHE)
  unset(CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO CACHE)
  unset(CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL CACHE)
  unset(CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO CACHE)
  unset(CMAKE_RC_FLAGS_MINSIZEREL CACHE)
  unset(CMAKE_RC_FLAGS_RELWITHDEBINFO CACHE)
endif()

# Default install dir:
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/.." CACHE STRING "" FORCE)

# compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Switches to enable examples and/or tests
option(EXAMPLES_ENABLED "enable examples" ON)
option(TESTS_ENABLED "enable tests" ON)

# Default install directory:
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/.." CACHE STRING "" FORCE)

# Internal include directory:
set(TOPE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "include path")
include_directories(${TOPE_INCLUDE_DIR})

# Scan for source and header files:
file(GLOB_RECURSE header_files include/*.h)
file(GLOB_RECURSE source_files src/*.c src/*.h)
source_group(header FILES ${header_files})
source_group(source FILES ${source_files})

# Define library:
add_library(tope ${header_files} ${source_files})
set_target_properties(tope PROPERTIES LINKER_LANGUAGE C)
if (UNIX)
  target_link_libraries(tope m)
endif()

# Install library and headers:
install(TARGETS tope ARCHIVE DESTINATION lib)
install(FILES ${header_files} DESTINATION include/tope)

# Tests:
if(TESTS_ENABLED)
  enable_testing()

  file(GLOB_RECURSE test_files tests/*.c)
  foreach(sname ${test_files})
    get_filename_component(name ${sname} NAME_WE)
    set (exename tope_${name})
    add_executable(${exename} ${sname})
    set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(${exename} tope)
    add_test(${exename} ${exename})
  endforeach()
endif()

# Examples:
if(EXAMPLES_ENABLED)
  file(GLOB_RECURSE example_files examples/*.c)
  foreach(sname ${example_files})
    get_filename_component(name ${sname} NAME_WE)
    set (exename tope_${name})
    add_executable(${exename} ${sname})
    set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(${exename} tope)
  endforeach()
endif()

# Executable:
add_executable(poly exe/main.c)
set_target_properties(poly PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(poly tope)
