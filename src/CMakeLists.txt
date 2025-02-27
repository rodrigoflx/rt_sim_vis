set(RT_SIM_VIS_SRC
    FluidSimulator.cpp
    main.cpp
    glad.c
)

PREPEND(RT_SIM_VIS_SRC)

set(RT_SIM_VIS_SRC ${RT_SIM_VIS_SRC}  PARENT_SCOPE)

add_library(glad STATIC glad.c)
add_executable(RTSimVisApp ${RT_SIM_VIS_SRC})

target_link_libraries(RTSimVisApp 
    glad
    glfw
    OpenGL::GL
)