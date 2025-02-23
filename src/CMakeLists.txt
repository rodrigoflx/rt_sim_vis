set(RT_SIM_VIS_SRC
    FluidSimulator.cpp
    main.cpp
)

PREPEND(RT_SIM_VIS_SRC)

set(RT_SIM_VIS_SRC ${RT_SIM_VIS_SRC}  PARENT_SCOPE)

add_executable(RTSimVisApp ${RT_SIM_VIS_SRC})
target_link_libraries(RTSimVisApp PRIVATE RTSimVis)