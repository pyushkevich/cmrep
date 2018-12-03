#ifndef __exp_approx__h_
#define __exp_approx__h_

/** Functions to compute an approximate exponential and its derivatives */
template <class TFloat>
inline TFloat 
exp_approx(TFloat x, TFloat f)
{
  TFloat z = 1.0 + (x * f) / 256.0;
  if(z < 0)
    return 0.0;

  z *= z; z *= z; z *= z; z *= z;
  z *= z; z *= z; z *= z; z *= z;

  return z;
}

template <class TFloat>
inline TFloat 
exp_approx(TFloat x, TFloat f, TFloat &d_fx)
{
  TFloat z = 1.0 + (x * f) / 256.0;
  if(z < 0)
    {
    d_fx = 0.0;
    return 0.0;
    }
    
  TFloat y = z;
  z *= z; z *= z; z *= z; z *= z;
  z *= z; z *= z; z *= z; z *= z;

  d_fx = f * z / y;

  return z;
}

template <class TFloat>
inline TFloat 
exp_approx(TFloat x, TFloat f, TFloat &d_fx, TFloat &d2_fx2)
{
  TFloat z = 1.0 + (x * f) / 256.0;
  TFloat scale = f / z;
  z *= z; z *= z; z *= z; z *= z;
  z *= z; z *= z; z *= z; z *= z;

  d_fx = z * scale;
  d2_fx2 = d_fx * scale * (255.0 / 256.0);

  return z;
}

#endif // __exp_approx__h_
