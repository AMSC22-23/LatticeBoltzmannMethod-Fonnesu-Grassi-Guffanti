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
add_tutorial(2_Pouseille.out "2_Pouseille.cpp")