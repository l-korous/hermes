project(ShockTubeProblem)

add_executable(${PROJECT_NAME} main.cpp ../euler_util.cpp numerical_flux.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
