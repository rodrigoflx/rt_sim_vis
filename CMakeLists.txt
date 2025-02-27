# Template based on the one presented in https://mjmorse.com/blog/cmake-template/

cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

project(RTSimVis VERSION 1.0)
set(RTSimVis VERSION 1.0)

cmake_policy(SET CMP0072 NEW)
set(OpenGL_GL_PREFERENCE GLVND)

include(CTest)
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})


find_package(glfw3 CONFIG REQUIRED)

# Function to prepend the subdirectory to source files in subdirectories
function(PREPEND var )
   set(listVar "")
   foreach(f ${${var}})
       list(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   endforeach(f)
   set(${var} "${listVar}" PARENT_SCOPE)
endfunction(PREPEND)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Put all required libraries here
find_package(Catch2 3 REQUIRED)


# Add subdirectories
add_subdirectory(src)
add_subdirectory(include)

# If testing, also add testing subdirectory
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

add_library(RTSimVis STATIC ${RT_SIM_VIS_SRC} ${RT_SIM_VIS_INC})

target_include_directories(RTSimVis PUBLIC /include)