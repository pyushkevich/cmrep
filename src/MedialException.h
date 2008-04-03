#ifndef __MedialException_h_
#define __MedialException_h_

#include <exception>
#include <stdexcept>

// Exception classes
class MedialModelException : public std::runtime_error 
{
public:
  MedialModelException(const char *text) : std::runtime_error(text) {}
};

class ModelIOException : public MedialModelException 
{
public:
  ModelIOException(const char *text) : MedialModelException(text) {}
};

#endif // __MedialException_h_
