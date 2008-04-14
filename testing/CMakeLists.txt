# Testing Commands
ENABLE_TESTING()

# Configure the parameter file
CONFIGURE_FILE(${CMREP_SOURCE_DIR}/testing/t001_param.txt ${CMREP_BINARY_DIR}/testing/t001_param.txt)

# Use the first test set
SET(TEST_TEMPLATE_BRUTE ${CMREP_SOURCE_DIR}/testing/t001_template_brute.cmrep)
SET(TEST_SUBJECT_BRUTE  ${CMREP_SOURCE_DIR}/testing/t001_subject_brute.cmrep)
SET(TEST_IMAGE_BINARY   ${CMREP_SOURCE_DIR}/testing/t001_img_binary.nii.gz)
SET(TEST_PARAM_FILE     ${CMREP_BINARY_DIR}/testing/t001_param.txt)

# Define the tests
ADD_TEST(TestBruteNoImage      ${CMREP_BINARY_DIR}/testpde DERIV1 ${TEST_SUBJECT_BRUTE})
ADD_TEST(TestBruteWithImage    ${CMREP_BINARY_DIR}/testpde DERIV2 ${TEST_SUBJECT_BRUTE} ${TEST_IMAGE_BINARY} ${TEST_PARAM_FILE})
