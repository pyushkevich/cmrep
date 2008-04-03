%module medialpde
%{
#include "ScriptInterface.h"
%}
%include "ScriptInterface.h"
%include exception.i

%exception {
  try {
    $function
  } catch(std::runtime_error &exc) {
    SWIG_exception(SWIG_SystemError, exc.what());
  } catch (...) {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}
