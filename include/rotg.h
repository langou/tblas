//
//  rotg.h
//
//  Purpose
//  =======
//
//  Generates rotation matrix elements (c,s) for a point (a,b) so that
//
//      (r,0) <- (c * a + s * b, -conj(s) * a + c * b)
//
//  where r = (x/|x|) sqrt(|x|^2 + |y|^2) if |x|!=0, or r = y if |x|=0.
//
//  Arguments
//  =========
//
//  n       length of vectors x and y
//
//  a       x-coordinate of point to rotate; on exit, contains r
//
//  b       y-coordinate of point to rotate
//
//  c       diagonal element of 2-by-2 rotation matrix (output only)
//
//  s       off-diagonal element of 2-by-2 rotation matrix (output only)
//

#ifndef __rotg__
#define __rotg__

#include <complex>
#include <cmath>

using std::complex;

namespace tblas
{
    template <typename T>
    void rotg(T &a, T &b, T &c, T &s)
    {
        using std::abs;
        using std::sqrt;
        using std::hypot;
        const T one(1.0);
        const T zero(0.0);
        T r=zero;
        T z=zero;
        T roe=(abs(a)>abs(b))?a:b;
        T scale=abs(a)+abs(b);
        if(scale==zero)
        {
            c=one;
            s=zero;
        }
        else
        {
            r=copysign(scale*hypot(a/scale,b/scale),roe);
            c=a/r;
            s=b/r;
            if(abs(a)>abs(b))
                z=s;
            else if(c!=zero)
                z=one/c;
            else
                z=one;
        }
        a=r;
        b=z;
    }

    template <typename T>
    void rotg(complex<T> &a,complex<T> b, T &c, complex<T> &s)
    {
        using std::abs;
        using std::sqrt;
        using std::norm;
        using std::hypot;
        const T zero(0.0);
        const complex<T> one(1.0);
        if(abs(a)==zero)
        {
            c=zero;
            s=one;
            a=b;
        }
        else
        {
            T scale=abs(a)+abs(b);
            T r=scale*hypot(abs(a/scale),abs(b/scale));
            complex<T> alpha=a/abs(a);
            c=abs(a)/r;
            s=alpha*conj(b)/r;
            a=alpha*r;
        }
    }
}
#endif
