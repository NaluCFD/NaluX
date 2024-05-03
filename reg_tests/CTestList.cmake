#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r)

# Standard performance test
function(add_test_p testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "performance")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_p)

# Regression test with single restart
function(add_test_r_rst testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE} && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.i -o ${testname}_rst.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}_rst ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_rst)

# Standard regression test with input data
function(add_test_r_t testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_t)

# Verification test with three resolutions
function(add_test_v3 testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R2.i -o ${testname}_R2.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v3)

# Verification test with two resolutions
function(add_test_v2 testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v2)

# Regression test that runs with different numbers of processes
function(add_test_r_np testname np)
    add_test(${testname}Np${np} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}Np${np}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}Np${np} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}Np${np}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname}Np${np} PROPERTIES TIMEOUT 1500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_np)

# Standard unit test
function(add_test_u testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${utest_ex_name}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "unit")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_u)


function(add_test_l testname np)
    add_test(${testname} sh -c "mpiexec -np ${np} --oversubscribe ${CMAKE_BINARY_DIR}/${nalu_ex_name} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/laboratory/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_CURRENT_SOURCE_DIR}/test_files/laboratory/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2500 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname}" LABELS "laboratory")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname}/mesh)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/laboratory/${testname}/mesh/${testname}.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/laboratory/${testname}/mesh)
endfunction(add_test_l)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(3dHex8ShockTube 2)

#=============================================================================
# Convergence tests
#=============================================================================
# nothing

#=============================================================================
# Unit tests
#=============================================================================
# nothing

#=============================================================================
# Performance tests
#=============================================================================
if(ENABLE_PERFORMANCE_TESTS)
# nothing
endif(ENABLE_PERFORMANCE_TESTS)

#=============================================================================
# Laboratory tests
#=============================================================================
if(ENABLE_LABORATORY_TESTS)
# nothing
endif(ENABLE_LABORATORY_TESTS)
