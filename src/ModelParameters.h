#ifndef THDMParameter_H
#define THDMParameter_H

#include <cmath>

#define Pi 3.1415926535897932384626433832795029
#define PI Pi
#define PI2 (PI*PI)
#define PIHalf (PI/2)

#define Alfa2 (Alfa*Alfa)

#define Sin(i) sin(i)
#define Csc(i) (1/sin(i))
#define Cos(i) cos(i)
#define Sec(i) (1/cos(i))
#define Tan(i) tan(i)
#define Cot(i) (1/tan(i))
#define Sqrt(i) (sqrt(i))
#define Power(a,b) (pow(a,b))
#define Complex(a,b) ComplexType(a,b)
#define Divergence (getdelta()) //For M2, we use MSbar scheme, so the whole expression will not be divergent, but it will explicitly contain symbol `Divergence`. In order to be consistent, we should use whatever LoopTools used for `Divergence`. By default it is 0.

#endif