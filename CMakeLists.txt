cmake_minimum_required(VERSION 3.5)

project(polytoop C)

if(UNIX)
  set(CMAKE_ASM_NASM_FLAGS "-f elf64 -DUNIX" CACHE STRING "" FORCE)

  set(CMAKE_C_FLAGS "-fvisibility=hidden -fPIC -no-pie -Wall -Wextra -std=c99 -pedantic" CACHE STRING "" FORCE)
  if(MEMCHECK)
    # valgind and AddressSanitizer don't mix
    set(CMAKE_C_FLAGS_DEBUG "-g ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
  else()
    # enable AddressSanitizer
    if (CMAKE_COMPILER_IS_GNUCC)
      # gcc
      set(CMAKE_C_FLAGS_DEBUG "-g -fsanitize=address -static-libasan ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
    else()
      # clang
      set(CMAKE_C_FLAGS_DEBUG "-g -fsanitize=address ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
    endif()
  endif()
  set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -s ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG -s ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
endif()

if (MSVC)
  set(CMAKE_ASM_NASM_FLAGS "-f win64 -DWIN64" CACHE STRING "" FORCE)

  # Replace dynamic runtime library (/MD) with static runtime library (/MT):
  set(variables CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE)
  foreach(variable ${variables})
    string(REGEX REPLACE "/MD" "/MT" var "${${variable}}")
    set(${variable} ${var} CACHE STRING "" FORCE)
  endforeach()

  # Suppress some warnings:
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4996 /wd6255 /wd6011")

  # No manifest:
  set(CMAKE_EXE_LINKER_FLAGS "/MANIFEST:NO" CACHE STRING "" FORCE)  # no EXE manifest
  set(CMAKE_SHARED_LINKER_FLAGS "/MANIFEST:NO" CACHE STRING "" FORCE)  # no DLL manifest
endif()

# Switches to enable examples and/or tests
option(EXAMPLES_ENABLED "enable examples" ON)
option(TESTS_ENABLED "enable tests" ON)

# Default install directory:
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/.." CACHE STRING "" FORCE)

# Internal include directory:
set(POLYTOOP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "include path")
include_directories(${POLYTOOP_INCLUDE_DIR})

# Scan for source and header files:
file(GLOB_RECURSE header_files include/*.h)
file(GLOB_RECURSE source_files src/*.c src/*.h)
source_group(header FILES ${header_files})
source_group(source FILES ${source_files})

# Define library:
add_library(polytoop ${header_files} ${source_files})
set_target_properties(polytoop PROPERTIES LINKER_LANGUAGE C)
if (UNIX)
  target_link_libraries(polytoop m)
endif()

# Install library and headers:
install(TARGETS polytoop ARCHIVE DESTINATION lib)
install(FILES ${header_files} DESTINATION include/polytoop)

# Tests:
if(TESTS_ENABLED)
  enable_testing()

  file(GLOB_RECURSE test_files tests/*.c)
  foreach(sname ${test_files})
    get_filename_component(name ${sname} NAME_WE)
    set (exename polytoop_${name})
    add_executable(${exename} ${sname})
    set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(${exename} polytoop)
    add_test(${exename} ${exename})
  endforeach()
endif()

# Examples:
if(EXAMPLES_ENABLED)
  file(GLOB_RECURSE example_files examples/*.c)
  foreach(sname ${example_files})
    get_filename_component(name ${sname} NAME_WE)
    set (exename polytoop_${name})
    add_executable(${exename} ${sname})
    set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(${exename} polytoop)
  endforeach()
endif()
