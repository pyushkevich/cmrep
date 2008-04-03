#ifndef __CodeTimer_h_
#define __CodeTimer_h_

#include <ctime>
#include <iostream>

class CodeTimer 
{
public:
  CodeTimer() 
    { tElapsed = 0.0; this->Start(); }

  void Start()
    { tStart = clock(); }
  
  void Stop()
    { tElapsed += (clock() - tStart) * 1.0 / CLOCKS_PER_SEC; }

  void Reset()
    { tElapsed = 0.0; }

  double Read()
    { return tElapsed; }

  double StopAndRead()
    { this->Stop(); return this->Read(); }

private:
  clock_t tStart;
  double tElapsed;
};


#endif // __CodeTimer_h_
