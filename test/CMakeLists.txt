macro(add_test_executable TARGET_NAME SOURCE_FILE)
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE LLALBM)
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/llalbm/include)
    message("${TARGET_NAME} uses include dir ${CMAKE_SOURCE_DIR}/llalbm/include")
endmacro()



# Adding reading executable
add_test_executable(LatticeReading.out "LatticeReadingTest.cpp")