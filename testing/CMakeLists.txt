# Configure the parameter file
CONFIGURE_FILE(${CMREP_SOURCE_DIR}/testing/t001_param.txt ${CMREP_BINARY_DIR}/testing/t001_param.txt)

# Use the first test set
SET(TEST_TEMPLATE_BRUTE ${CMREP_SOURCE_DIR}/testing/t001_template_brute.cmrep)
SET(TEST_SUBJECT_BRUTE  ${CMREP_SOURCE_DIR}/testing/t001_subject_brute.cmrep)
SET(TEST_SUBJECT_PDE    ${CMREP_SOURCE_DIR}/testing/caudate_pde_model.cmrep)
SET(TEST_IMAGE_BINARY   ${CMREP_SOURCE_DIR}/testing/t001_img_binary.nii.gz)
SET(TEST_IMAGE_CAUDATE  ${CMREP_SOURCE_DIR}/testing/caudate_target.nii.gz)
SET(TEST_PARAM_FILE     ${CMREP_BINARY_DIR}/testing/t001_param.txt)

# Define the tests
ADD_TEST(TestSparseSolver      ${CMREP_BINARY_DIR}/cmrep_test SPARSE)
ADD_TEST(TestBruteNoImage      ${CMREP_BINARY_DIR}/cmrep_test DERIV1 ${TEST_SUBJECT_BRUTE})
ADD_TEST(TestBruteWithImage    ${CMREP_BINARY_DIR}/cmrep_test DERIV2 ${TEST_SUBJECT_BRUTE} ${TEST_IMAGE_BINARY} ${TEST_PARAM_FILE})
ADD_TEST(TestPDENoImage        ${CMREP_BINARY_DIR}/cmrep_test DERIV1 ${TEST_SUBJECT_PDE})
ADD_TEST(TestPDEWithImage      ${CMREP_BINARY_DIR}/cmrep_test DERIV2 ${TEST_SUBJECT_PDE} ${TEST_IMAGE_CAUDATE})

# Geodesic shooting tests
SET(TEST_GSHOOT_DIR ${CMREP_SOURCE_DIR}/testing/shooting)
ADD_TEST(ShootingRegression    ${CMREP_BINARY_DIR}/shooting_test ${TEST_GSHOOT_DIR}/shoot_regression.mat)
ADD_TEST(ShootingA             ${CMREP_BINARY_DIR}/lmshoot 
                                 -m ${TEST_GSHOOT_DIR}/shooting_test_2d_template.vtk 
                                    ${TEST_GSHOOT_DIR}/shooting_test_2d_target.vtk 
                                 -o ${CMREP_BINARY_DIR}/res.vtk 
                                 -s 0.2 -l 20 -d 2 -a A -i 20)
ADD_TEST(ShootingG             ${CMREP_BINARY_DIR}/lmshoot 
                                 -m ${TEST_GSHOOT_DIR}/shooting_test_2d_template.vtk 
                                    ${TEST_GSHOOT_DIR}/shooting_test_2d_target.vtk 
                                 -o ${CMREP_BINARY_DIR}/res.vtk 
                                 -s 0.2 -l 20 -d 2 -a G -i 20)
