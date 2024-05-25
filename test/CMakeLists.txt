macro(add_test_executable TARGET_NAME SOURCE_FILE)
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE LLALBM)
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/llalbm/include)
    message("${TARGET_NAME} uses include dir ${CMAKE_SOURCE_DIR}/llalbm/include")
endmacro()



# Adding reading executable
add_test_executable(LatticeReading.out "unit/LatticeReadingTest.cpp")
add_test_executable(Lid.out "unit/LiddrivenTest.cpp")
add_test_executable(InletUpdate.out "unit/InletUpdates.cpp")
add_test_executable(LidExec.out "unit/LidExecTest.cpp")
add_test_executable(GenerationTest.out "unit/GenerationTest.cpp")
add_test_executable(LidFromGeneration.out "unit/LidFromGeneration.cpp")
