macro(add_tutorial TARGET_NAME SOURCE_FILE)
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE LLALBM)
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/llalbm/include)
    message("${TARGET_NAME} uses include dir ${CMAKE_SOURCE_DIR}/llalbm/include")
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
    endif()
endmacro()



# Adding reading executable
add_tutorial(1_LidDrivenCavity.out "1_LidDrivenCavity.cpp")
add_tutorial(2_Parabolic.out "2_Parabolic.cpp")
add_tutorial(3_Vortex.out "3_Vortex.cpp")
add_tutorial(4_ComplexGeometry.out "4_ComplexGeometry.cpp")
add_tutorial(5_RandomGeneration.out "5_RandomGeneration.cpp")
add_tutorial(6_FlowAnalysis.out "6_FlowAnalysis.cpp")
add_tutorial(7_FlatBottomFoil.out "7_FlatBottomFoil.cpp")
add_tutorial(8_1_NACAProfileFlight.out "8_1_NACAProfileFlight.cpp")
add_tutorial(8_2_NACAProfileTakeOff.out "8_2_NACAProfileTakeOff.cpp")
add_tutorial(9_Conduct.out "9_Conduct.cpp")