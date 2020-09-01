//
//  rotm.h
//
//  Purpose
//  =======
//
//  Applies a modified Givens rotation matrix H to a set of points (x,y).
//
//  Arguments
//  =========
//
//  n       length of vectors x and y
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  y       vector of length n
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//
//  H       rotation matrix from rotmg stored as array of length 4
//
//  flag    mode of rotation matrix provided on return from rotmg
//

#ifndef __rotm__
#define __rotm__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void rotm(size_t n, T *x, ptrdiff_t incx, T *y, ptrdiff_t incy, T *H, int flag)
    {
        if((incx==incy)&&(incx==1))
        {
            if(flag==-1)
            {
                for(size_t i=0;i<n;i++)
                {
                    T w=x[i];
                    T z=y[i];
                    x[i]=w*H[0]+z*H[2];
                    y[i]=w*H[1]+z*H[3];
                }
            }
            else if(flag==0)
            {
                for(size_t i=0;i<n;i++)
                {
                    T w=x[i];
                    T z=y[i];
                    x[i]=w+z*H[2];
                    y[i]=w*H[1]+z;
                }
            }
            else if(flag==1)
            {
                for(size_t i=0;i<n;i++)
                {
                    T w=x[i];
                    T z=y[i];
                    x[i]=w*H[0]+z;
                    y[i]=-w+z*H[3];
                }
            }
        }
        else
        {
            size_t ix=(incx>0)?0:(1-n)*incx;
            size_t iy=(incy>0)?0:(1-n)*incy;
            if(flag==-1)
            {
                for(size_t i=0;i<n;i++)
                {
                    T w=x[ix];
                    T z=y[iy];
                    x[ix]=w*H[0]+z*H[2];
                    y[iy]=w*H[1]+z*H[3];
                    ix+=incx;
                    iy+=incy;
                }
            }
            else if(flag==0)
            {
                for(size_t i=0;i<n;i++)
                {
                    T w=x[ix];
                    T z=y[iy];
                    x[ix]=w+z*H[2];
                    y[iy]=w*H[1]+z;
                    ix+=incx;
                    iy+=incy;
                }
            }
            else if(flag==1)
            {
                for(size_t i=0;i<n;i++)
                {
                    T w=x[ix];
                    T z=y[iy];
                    x[ix]=w*H[0]+z;
                    y[iy]=-w+z*H[3];
                    ix+=incx;
                    iy+=incy;
                }
            }
            
        }
    }
}
#endif
