//
//  rotmg.h
//
//  Purpose
//  =======
//
//  Constructs the modified Givens transformation matrix H.
//
//  The Givens rotation matrix H is constructed that zeros the second
//  component of the vector (sqrt(d1)*x1,sqrt(d2)*y1).
//
//  The value of flag indicats the form of H as follows:
//
//          flag=-1          flag=0          flag=1         flag=-2
//        ------------    ------------    -------------    ---------
//        [ h11  h12 ]    [  1   h12 ]    [ h11    1  ]    [ 1   0 ]
//        [ h21  h22 ]    [ h21   1  ]    [ -1    h22 ]    [ 0   1 ]
//
//  Returns
//  =======
//
//  flag indicating mode of rotation matrix
//
//  Arguments
//  =========
//
//  d1      scalar (see above)
//
//  d2      scalar (see above)
//
//  x1      scalar (see above)
//
//  x2      scalar (see above)
//
//  H       rotation matrix stored as array of length 4
//

#ifndef __rotmg__
#define __rotmg__

#include <cmath>

namespace tblas
{
    template <typename T>
    int rotmg(T &d1, T &d2, T &x1, T &y1, T *H)
    {
        using std::abs;
        const T zero(0.0);
        const T one(1.0);
        const T gam(4096.0);
        const T gamsq(gam*gam);
        const T rgamsq(one/gamsq);
        T r=zero;
        T u=zero;
        T p1=zero;
        T p2=zero;
        T q1=zero;
        T q2=zero;
        T h11=zero;
        T h12=zero;
        T h21=zero;
        T h22=zero;
        int flag=0;
        T stemp=zero;
        
        if(d1<zero)
        {
            flag=-1;
            h11=zero;
            h12=zero;
            h21=zero;
            h22=zero;
            d1=zero;
            d2=zero;
            x1=zero;
        }
        else
        {
            p2=d2*y1;
            if(p2==zero)
            {
                flag=-2;
                return flag;
            }
            p1=d1*x1;
            q2=p2*y1;
            q1=p1*x1;
            if(abs(q1)>abs(q2))
            {
                h21=-(y1)/x1;
                h12=p2/p1;
                u=one-h12*h21;
                if(u>zero)
                {
                    flag=0;
                    d1/=u;
                    d2/=u;
                    x1*=u;
                }
            }
            else
            {
                if(q2<zero)
                {
                    flag=-1;
                    h11=zero;
                    h12=zero;
                    h21=zero;
                    h22=zero;
                    d1=zero;
                    d2=zero;
                    x1=zero;
                }
                else
                {
                    flag=1;
                    h11=p1/p2;
                    h22=x1/y1;
                    u=one+h11*h22;
                    stemp=d2/u;
                    d2=d1/u;
                    d1=stemp;
                    x1=y1*u;
                }
            }
            if(d1!=zero)
            {
                while((d1<=rgamsq)||(d1>=gamsq))
                {
                    if(flag==0)
                    {
                        h11=one;
                        h22=one;
                        flag=-1;
                    }
                    else
                    {
                        h21=-one;
                        h12=one;
                        flag=-1;
                    }
                    if(d1<=rgamsq)
                    {
                        r=gam;
                        d1*=r*r;
                        x1/=gam;
                        h11/=gam;
                        h12/=gam;
                    }
                    else
                    {
                        r=gam;
                        d1/=r*r;
                        x1*=gam;
                        h11*=gam;
                        h12*=gam;
                    }
                }
            }
            if(d2!=zero)
            {
                while((abs(d2)<=rgamsq)||(abs(d2)>=gamsq))
                {
                    if(flag==0)
                    {
                        h11=one;
                        h22=one;
                        flag=-1;
                    }
                    else
                    {
                        h21=-one;
                        h12=one;
                        flag=-1;
                    }
                    if(abs(d2)<=rgamsq)
                    {
                        r=gam;
                        d2*=r*r;
                        h21/=gam;
                        h22/=gam;
                    }
                    else
                    {
                        r=gam;
                        d2/=r*r;
                        h21*=gam;
                        h22*=gam;
                    }
                }
            }
        }
        if(flag<0)
        {
            H[0]=h11;
            H[1]=h21;
            H[2]=h12;
            H[3]=h22;
        }
        else if(flag==0)
        {
            H[1]=h21;
            H[2]=h12;
        }
        else
        {
            H[0]=h11;
            H[3]=h22;
        }
        return flag;
    }
}
#endif
