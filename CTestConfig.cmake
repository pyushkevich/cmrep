## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "cmrep")
set(CTEST_NIGHTLY_START_TIME "04:05:00 EST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "itksnap.org")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=cmrep")
set(CTEST_DROP_SITE_CDASH TRUE)
