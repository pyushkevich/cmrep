// Setup printing of stack trace on segmentation faults. This only
// works on select GNU systems
#if defined(__GNUC__) && !defined(__CYGWIN__) && !defined(__APPLE__) && !defined(sun) && !defined(WIN32)

#include <signal.h>
#include <execinfo.h>
#include <iostream>

using namespace std;

void SegmentationFaultHandler(int sig)
{
  cerr << "*************************************" << endl;
  cerr << "ITK-SNAP: Segmentation Fault!   " << endl;
  cerr << "BACKTRACE: " << endl;
  void *array[50];
  int nsize = backtrace(array, 50);
  backtrace_symbols_fd(array, nsize, 2);
  cerr << "*************************************" << endl;
  exit(-1);
}

void SetupSignalHandlers()
{
  signal(SIGSEGV, SegmentationFaultHandler);
}

#else

void SetupSignalHandlers()
{
  // Nothing to do!
}

#endif



