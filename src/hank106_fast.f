c
c
c
        subroutine hank106init(rk7,rmin,rmax,w,keep,ninterv)
        implicit real *8 (a-h,o-z)
        complex *16 rk,h0,h1,ima,z,rk7,u0,u07,com
        dimension w(*),rea(2),rws(10)
        integer iws(20)
        
c
        equivalence (rea(1),com),(rws(1),iws(1))
c
        data ima/(0.0d0,1.0d0)/
c
c        This subroutine evaluates the Hankel functions 
c        H_0^1, H^1_1 of a complex argument, the argument 
c        living on a ray. The subroutine you are looking at is
c        the initialization subroutine; the evaluation 
c        subroutine is hank106 (see).
c
c
c         8/28/02 - added ninterv to calling sequence  (LG)
c
c
c   PLEASE NOTE THAT THE USE OF THIS SUBROUTIONE IS NOT COMPLETELY
C   STRAIGHTFORWARD: A STRAIGHTFORWARD SUBROUTINE TO USE IS HANK103
C   (SEE). HOWEVER, THIS SUBROUTINE IS ABOUT 4 TIMES FASTER THAN 
C   HANK103.
c
c        Recommended pairs rmin,rmax (assuming that rk \sim 1):
c
c
c
c        rmin       rmax
c
c        0.062      0.125
c        0.031      0.062
c        0.062       0.125
c        0.125       0.25
c        0.25        0.5
c        0.5         1
c        1           2
c        2           5
c        5           10
c        10          100
c        10          200
c
c  GENERALLY, WITH ABS(RMIN*RK) > 10, THERE IS NO SHARP LIMIT ON
C  RMAX. HOWEVER, FOR SUFFICIENTLY LARGE (RMAX-RMIN)*RK, THE 
C  code loses speed due to caching problems. THE DETERIORATION
C  BECOMES NOTICEABLE AT SOME POINT AFTER |(RMAX-RMIN)*RK| > 100
C  (ON THE PENTIUM-IV DESKTOP).
C
c
c                   Input parameters:
c
c  rk7 - the Helmholtz coefficient
c  rmin - the minimum r for which this subroutine will evaluate 
c       the Hankel functions
c  rmax - the minimum r for which this subroutine will evaluate 
c       the Hankel functions
c
c                   Output parameters:
c
c  w - contains various data to be used by the entry hank106 (see)
c  keep - the first keep elements of the array w should not be 
c       changed between the call to this entry, and the subsequent
c       calls to the entry hank106.
c  ninterv - number of intervals used in equispaced
c       subdivision of current interval
c
c
        rk=rk7
c
        dfool=rk/abs(rk)
        if(abs(dfool) .lt. 1.0d-30) rk=rk+2*1.0d-30*abs(rk)
c
        d=rk
c
        x1=d*rmin
        x2=d*rmax
c
c        initialize the evaluation of H_0, H_1
c
        n=11
        ddd=(abs(rmax*rk)-abs(rmin*rk))*2
        i=ddd
c
 
ccc        call prin2('ddd=*',ddd,1)
        if(i .lt. 10) i=10
        ninterv=i
c
ccc        call prinf('ninterv as calculated*',ninterv,1)
c
        d=rk
        d2=-ima*rk
        coef=d2/d
c
c       allocate memory for the initialization
c
        icenters=21
        lcenters=ninterv*2+4
c
        ih0s=icenters+lcenters
        lh0s=ninterv+4
        lh0s=lh0s*2
c
        ih1s=ih0s+lh0s
        lh1s=ninterv+4
        lh1s=lh1s*2
c
        ih0derss=ih1s+lh1s
        lh0derss=(ninterv*n+4)*2
c
        ih1derss=ih0derss+lh0derss
        lh1derss=(ninterv*n+4)*2
c
        keep=ih1derss+lh1derss
c
        call hank106ini(coef,x1,x2,ninterv,n,
     1      w(icenters),w(ih0derss),w(ih0s),w(ih1s),h,
     2      w(ih1derss),u07)
c
        u0=1/u07
c
c       store in the beginning of the array w various types of data 
c
        ix1=1
        ih=2
        iu0=3
        in=5
c
        irk=6
c
        w(ix1)=x1
        w(ih)=h
        w(ih)=1/h
c
        com=u0
        w(iu0)=rea(1)
        w(iu0+1)=rea(2)
c
        w(in)=n+0.1
c
        com=rk
c
        w(irk)=rea(1)
        w(irk+1)=rea(2)
c
c        store integer data in the array w
c
        iws(1)=ih0derss
        iws(2)=ih0s
        iws(3)=ih1derss
        iws(4)=ih1s

c
        do 3200 j=1,8
c
        w(10+j)=rws(j)
 3200 continue
c
        return
        end
c
c
c
c
c
        subroutine hank106ini(coef,x1,x2,ninterv,n,
     1      centers,h0derss,h0s,h1s,h,h1derss,u0)
        implicit real *8 (a-h,o-z)
        complex *16 h0s(1),h0derss(n,1),h1derss(n,1),ima,u0,
     1      us(22),h1s(1)
        dimension centers(2,1)
c
        data ima/(0.0d0,1.0d0)/
c 
c        construct the subintervals
c
        h=(x2-x1)/ninterv
c
        do 1200 i=1,ninterv
c
        ab1i=(i-1)*h+x1
        ab2i=(i-1)*h+x1 +h
        centers(1,i)=(ab2i+ab1i)/2
        centers(2,i)=coef*centers(1,i)
 1200 continue
c
c        construct the values of Hankel functions and their
c        derivatibes at the centers
c
        do 1400 i=1,ninterv
c
        call hank0ders(centers(1,i),n,h0s(i),h1s(i),
     1      h0derss(1,i),h1derss(1,i) )
 1400 continue
c
c       scale them things by factorials and by complex powers
c
        u0=1+ima*coef
        u0=u0/abs(u0)
c
        us(1)=u0
        do 1500 i=1,20
c
        us(i+1)=us(i)*u0
 1500 continue
c
        do 1800 i=1,ninterv
        fact=1
        do 1600 j=1,n-1
        h0derss(j,i)=h0derss(j,i)*fact * us(j)
        h1derss(j,i)=h1derss(j,i)*fact * us(j)
c
        fact=fact/(j+1)
 1600 continue
c
 1800  continue
c
        return
        end
c
c
c
c
c
        subroutine hank0ders(z,n,h0,h1,h0ders,h1ders)
        implicit real *8 (a-h,o-z)
        complex *16 z,h0,h0ders(1),h1,h1ders(1)
c
        data ifexpon/1/
c 
c        evaluate h0 and h1
c
        call hank103(z,h0,h1,ifexpon)
c
        h0ders(1)=-h1
        h0ders(2)=-(h0ders(1)/z+h0)
        h0ders(3)=-(2*h0ders(2)+h0ders(1)*z+h0)/z
        h0ders(4)=-(3*h0ders(3)+h0ders(2)*z+2*h0ders(1))/z
c
        if(n .le. 4) return
c
        do 1400 m=2,n-2
c
        h0ders(m+2)=-( (m+1)*h0ders(m+1)+z*h0ders(m)+
     1      m*h0ders(m-1) )/z
 1400 continue
c
        do 1600 i=1,n-1
c
        h1ders(i)=-h0ders(i+1)
 1600 continue
c
        return
        end
c
c
c
c
c

        subroutine hank106datagen_r(w,lw,ier)
        implicit real *8 (a-h,o-z)
        integer istart(29),nintervec(28)
        dimension w(lw),ab(2,28)
        complex *16 rk,h0,h1,z,rksav
        data nab/28/
        data rmin/1.0d-6/
        data rmax/200/
c
c     INPUT PARAMETERS:    ----> Now hidden as local vars.....
c
c       create top level (dyadic) intervals for hank106init, which 
c       then uses equisized subintervals to precompute interpolation
c       polynomials
c
c     rk (complex *16) frequency parameter
c
c     rmin, rmax (real *8)  desired range of argument to hank106
c                           [rmin*rk,...,rmax*rk]
c 
c     ab(2,nab) (real *8)   blank array of length 2*nab
c     w(lw)     (real *8)   work array of length lw
c
c     OUTPUT PARAMETERS:
c
c     ninterval (integer *4)   number of subintervals created
c     ab(2,ninterval) (real *8) boundary of ith interval is
c                               (ab(1,i),ab(2,i))
c     nintervec(ninterval)    nuomber of equispaced subintervals
c                             used for ith interval
c     istart  (integer *4)     istart(i) is pointer into workspace for
c                              data pertaining to ith interval
c     ier (integer *4) error flag
c             ier = 0 upon normal execution.
c             ier = 1 if length (nab) of array ab is of 
c                     insufficient length 
c             ier = 2 if length (lw) of workspace w is of 
c                     insufficient length 
c-----------------------------------------------------------------
c
c
        ier = 0
        rk = 1.0d0
        ninterval = 1
        istart(1) = 1
        rmaxloc = rmin
        rksav = rk/cdabs(rk)
        rminsav = rmin
        do i = 1,28
           rminloc = rmaxloc 
           rmaxloc = 2*rminloc 
           if (rminloc.gt.100) rmaxloc = rminloc+100
           ab(1,i) = rminloc
           ab(2,i) = rmaxloc
           call hank106init(rksav,rminloc,rmaxloc,w(istart(i)),
     1        keep,ninterv)
           nintervec(i) = ninterv
           istart(i+1) = istart(i) + keep + 1
           if (rmaxloc.ge.rmax) goto 1111
           if (istart(i+1).gt.lw) then
              ier = 2
              return
           endif
           ninterval = ninterval + 1
        enddo
1111    continue

        return
        end
c
c
c
c
c
        subroutine hank106_r_h0(r,h0,ifexpon,intnum,w,lw)
        implicit real *8 (a-h,o-z)
        real *8 r
        complex *16 h0,z,h1
        integer ifexpon,ninterval,nintervec(28)
        integer istart(28),ijws(4,28)
        real *8 w(lw)
        real *8 ab(2,28)
        data ninterval/28/
        data nintervec/10,10,10,10,10,10,10,10,10,10,
     1   10,10,10,10,10,10,10,10,10,10,
     2   10,10,10,16,33,67,134,200/
        data ijws/100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            100,44,328,72,
     4            100,44,328,72,
     5            100,44,328,72,
     6            100,44,328,72,
     7            100,44,328,72,
     8            100,44,328,72,
     9            100,44,328,72,
     *            100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            100,44,328,72,
     4            100,44,328,72,
     5            100,44,328,72,
     6            100,44,328,72,
     7            100,44,328,72,
     8            100,44,328,72,
     9            100,44,328,72,
     *            100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            136,56,496,96,
     4            238,90,972,164,
     5            442,158,1924,300,
     6            844,292,3800,568,
     7            1240,424,5648,832/
        data ab/
     1     0.1000000000000000E-05,  0.2000000000000000E-05,
     1     0.2000000000000000E-05,  0.4000000000000000E-05,
     1     0.4000000000000000E-05,  0.8000000000000000E-05,
     1     0.8000000000000000E-05,  0.1600000000000000E-04,
     1     0.1600000000000000E-04,  0.3200000000000000E-04,
     1     0.3200000000000000E-04,  0.6400000000000000E-04,
     1     0.6400000000000000E-04,  0.1280000000000000E-03,
     1     0.1280000000000000E-03,  0.2560000000000000E-03,
     1     0.2560000000000000E-03,  0.5120000000000000E-03,
     1     0.5120000000000000E-03,  0.1024000000000000E-02,
     1     0.1024000000000000E-02,  0.2048000000000000E-02,
     1     0.2048000000000000E-02,  0.4096000000000000E-02,
     1     0.4096000000000000E-02,  0.8192000000000000E-02,
     1     0.8192000000000000E-02,  0.1638400000000000E-01,
     1     0.1638400000000000E-01,  0.3276800000000000E-01,
     1     0.3276800000000000E-01,  0.6553600000000000E-01,
     1     0.6553600000000000E-01,  0.1310720000000000E+00,
     1     0.1310720000000000E+00,  0.2621440000000000E+00,
     1     0.2621440000000000E+00,  0.5242880000000000E+00,
     1     0.5242880000000000E+00,  0.1048576000000000E+01,
     1     0.1048576000000000E+01,  0.2097152000000000E+01,
     1     0.2097152000000000E+01,  0.4194304000000000E+01,
     1     0.4194304000000000E+01,  0.8388608000000000E+01,
     1     0.8388608000000000E+01,  0.1677721600000000E+02,
     1     0.1677721600000000E+02,  0.3355443200000000E+02,
     1     0.3355443200000000E+02,  0.6710886400000000E+02,
     1     0.6710886400000000E+02,  0.1342177280000000E+03,
     1     0.1342177280000000E+03,  0.2342177280000000E+03/
       data istart/1,559,1117,1675,2233,2791,3349,3907,
     1    4465,5023,5581,6139,6697,7255,7813,8371,8929,
     2    9487,10045,10603,11161,11719,12277,12835,13693,15401,18809,
     3    25567/

c
c
c
c-------------------------------------------------------------
c
c       determine subinterval and call hank106.
c
        call findinte(r,ab,ninterval,intnum)
        i = intnum
        i1 = istart(i)
        if (intnum.le.ninterval.and.intnum.ge.1) then
          call hank106eva(r,w(i1),w(i1+ijws(1,i)),
     1     w(i1+ijws(2,i)),h0,w(i1+1),nintervec(i))
        else
           z = r
           call hank103(z,h0,h1,ifexpon)
        endif
        return
        end
c
c
c
c
c
c
c
c
c
c
c
        subroutine hank106_r_h1(r,h1,ifexpon,intnum,w,lw)
        implicit real *8 (a-h,o-z)
        real *8 r
        complex *16 h0,z,h1
        integer ifexpon,ninterval,nintervec(28)
        integer istart(28),ijws(4,28)
        real *8 w(lw)
        real *8 ab(2,28)
        data ninterval/28/
        data nintervec/10,10,10,10,10,10,10,10,10,10,
     1   10,10,10,10,10,10,10,10,10,10,
     2   10,10,10,16,33,67,134,200/
        data ijws/100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            100,44,328,72,
     4            100,44,328,72,
     5            100,44,328,72,
     6            100,44,328,72,
     7            100,44,328,72,
     8            100,44,328,72,
     9            100,44,328,72,
     *            100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            100,44,328,72,
     4            100,44,328,72,
     5            100,44,328,72,
     6            100,44,328,72,
     7            100,44,328,72,
     8            100,44,328,72,
     9            100,44,328,72,
     *            100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            136,56,496,96,
     4            238,90,972,164,
     5            442,158,1924,300,
     6            844,292,3800,568,
     7            1240,424,5648,832/
        data ab/
     1     0.1000000000000000E-05,  0.2000000000000000E-05,
     1     0.2000000000000000E-05,  0.4000000000000000E-05,
     1     0.4000000000000000E-05,  0.8000000000000000E-05,
     1     0.8000000000000000E-05,  0.1600000000000000E-04,
     1     0.1600000000000000E-04,  0.3200000000000000E-04,
     1     0.3200000000000000E-04,  0.6400000000000000E-04,
     1     0.6400000000000000E-04,  0.1280000000000000E-03,
     1     0.1280000000000000E-03,  0.2560000000000000E-03,
     1     0.2560000000000000E-03,  0.5120000000000000E-03,
     1     0.5120000000000000E-03,  0.1024000000000000E-02,
     1     0.1024000000000000E-02,  0.2048000000000000E-02,
     1     0.2048000000000000E-02,  0.4096000000000000E-02,
     1     0.4096000000000000E-02,  0.8192000000000000E-02,
     1     0.8192000000000000E-02,  0.1638400000000000E-01,
     1     0.1638400000000000E-01,  0.3276800000000000E-01,
     1     0.3276800000000000E-01,  0.6553600000000000E-01,
     1     0.6553600000000000E-01,  0.1310720000000000E+00,
     1     0.1310720000000000E+00,  0.2621440000000000E+00,
     1     0.2621440000000000E+00,  0.5242880000000000E+00,
     1     0.5242880000000000E+00,  0.1048576000000000E+01,
     1     0.1048576000000000E+01,  0.2097152000000000E+01,
     1     0.2097152000000000E+01,  0.4194304000000000E+01,
     1     0.4194304000000000E+01,  0.8388608000000000E+01,
     1     0.8388608000000000E+01,  0.1677721600000000E+02,
     1     0.1677721600000000E+02,  0.3355443200000000E+02,
     1     0.3355443200000000E+02,  0.6710886400000000E+02,
     1     0.6710886400000000E+02,  0.1342177280000000E+03,
     1     0.1342177280000000E+03,  0.2342177280000000E+03/
       data istart/1,559,1117,1675,2233,2791,3349,3907,
     1    4465,5023,5581,6139,6697,7255,7813,8371,8929,
     2    9487,10045,10603,11161,11719,12277,12835,13693,15401,18809,
     3    25567/

c
c
c
c-------------------------------------------------------------
c
c       determine subinterval and call hank106.
c
        call findinte(r,ab,ninterval,intnum)
        i = intnum
        i1 = istart(i)
        if (intnum.le.ninterval.and.intnum.ge.1) then
          call hank106eva(r,w(i1),w(i1+ijws(3,i)),
     1     w(i1+ijws(4,i)),h1,w(i1+1),nintervec(i))
        else
           z = r
           call hank103(z,h0,h1,ifexpon)
        endif
        return
        end
c
c
c
c
c
c
c
        subroutine hank106_r_h01(r,h0,h1,ifexpon,intnum,w,lw)
        implicit real *8 (a-h,o-z)
        real *8 r
        complex *16 h0,z,h1
        integer ifexpon,ninterval,nintervec(28)
        integer istart(28),ijws(4,28)
        real *8 w(lw)
        real *8 ab(2,28)
        data ninterval/28/
        data nintervec/10,10,10,10,10,10,10,10,10,10,
     1   10,10,10,10,10,10,10,10,10,10,
     2   10,10,10,16,33,67,134,200/
        data ijws/100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            100,44,328,72,
     4            100,44,328,72,
     5            100,44,328,72,
     6            100,44,328,72,
     7            100,44,328,72,
     8            100,44,328,72,
     9            100,44,328,72,
     *            100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            100,44,328,72,
     4            100,44,328,72,
     5            100,44,328,72,
     6            100,44,328,72,
     7            100,44,328,72,
     8            100,44,328,72,
     9            100,44,328,72,
     *            100,44,328,72,
     1            100,44,328,72,
     2            100,44,328,72,
     3            136,56,496,96,
     4            238,90,972,164,
     5            442,158,1924,300,
     6            844,292,3800,568,
     7            1240,424,5648,832/
        data ab/
     1     0.1000000000000000E-05,  0.2000000000000000E-05,
     1     0.2000000000000000E-05,  0.4000000000000000E-05,
     1     0.4000000000000000E-05,  0.8000000000000000E-05,
     1     0.8000000000000000E-05,  0.1600000000000000E-04,
     1     0.1600000000000000E-04,  0.3200000000000000E-04,
     1     0.3200000000000000E-04,  0.6400000000000000E-04,
     1     0.6400000000000000E-04,  0.1280000000000000E-03,
     1     0.1280000000000000E-03,  0.2560000000000000E-03,
     1     0.2560000000000000E-03,  0.5120000000000000E-03,
     1     0.5120000000000000E-03,  0.1024000000000000E-02,
     1     0.1024000000000000E-02,  0.2048000000000000E-02,
     1     0.2048000000000000E-02,  0.4096000000000000E-02,
     1     0.4096000000000000E-02,  0.8192000000000000E-02,
     1     0.8192000000000000E-02,  0.1638400000000000E-01,
     1     0.1638400000000000E-01,  0.3276800000000000E-01,
     1     0.3276800000000000E-01,  0.6553600000000000E-01,
     1     0.6553600000000000E-01,  0.1310720000000000E+00,
     1     0.1310720000000000E+00,  0.2621440000000000E+00,
     1     0.2621440000000000E+00,  0.5242880000000000E+00,
     1     0.5242880000000000E+00,  0.1048576000000000E+01,
     1     0.1048576000000000E+01,  0.2097152000000000E+01,
     1     0.2097152000000000E+01,  0.4194304000000000E+01,
     1     0.4194304000000000E+01,  0.8388608000000000E+01,
     1     0.8388608000000000E+01,  0.1677721600000000E+02,
     1     0.1677721600000000E+02,  0.3355443200000000E+02,
     1     0.3355443200000000E+02,  0.6710886400000000E+02,
     1     0.6710886400000000E+02,  0.1342177280000000E+03,
     1     0.1342177280000000E+03,  0.2342177280000000E+03/
       data istart/1,559,1117,1675,2233,2791,3349,3907,
     1    4465,5023,5581,6139,6697,7255,7813,8371,8929,
     2    9487,10045,10603,11161,11719,12277,12835,13693,15401,18809,
     3    25567/

c
c
c
c-------------------------------------------------------------
c
c       determine subinterval and call hank106.
c
        call findinte(r,ab,ninterval,intnum)
        i = intnum
        i1 = istart(i)
        if (intnum.le.ninterval.and.intnum.ge.1) then
          call hank106eva2(r,w(i1),w(i1+ijws(1,i)),
     1     w(i1+ijws(2,i)),w(i1+ijws(3,i)),
     1     w(i1+ijws(4,i)),h0,h1,w(i1+1),nintervec(i))
        else
           z = r
           call hank103(z,h0,h1,ifexpon)
        endif
        return
        end
c
c
c
c
c
        subroutine hank106eva(r,x1,h0derss,
     1      h0s,h0,h,ninterv)
        implicit real *8 (a-h,o-z)
        real *8 h0s(2,ninterv),h0derss(2,11,ninterv),h0(2)
        real *8 x1,h
c
c 
c   input:
c   ninterv - number of intervals used in equispaced
c       subdivision of current interval
c
c-----------------------------------
c
c       find the subinterval where the point z lives
c
        i=(r-x1)*h +1
        if (i.lt.0) then 
            i = 1
        else if (i.gt.ninterv) then 
            i = ninterv
        endif

        t = (r-x1) - (2*i-1)/h/2
c
c        evaluate the functions h0 and h1 at the point z
c
c
        h0(:)=(((((((((h0derss(:,10,i)*t+
     2      h0derss(:,9,i))*t+h0derss(:,8,i) ) 
     1    * t+h0derss(:,7,i))*t+h0derss(:,6,i))*t+h0derss(:,5,i))*t 
     2    +h0derss(:,4,i))*t+h0derss(:,3,i))*t+h0derss(:,2,i))
     3    *t+h0derss(:,1,i) ) * t + h0s(:,i)
c
c
        return
        end
c
c
c
c
c
c
c
        subroutine hank106eva2(r,x1,h0derss,
     1      h0s,h1derss,h1s,h0,h1,h,ninterv)
        implicit real *8 (a-h,o-z)
        real *8 h0s(2,ninterv),h0derss(2,11,ninterv),h0(2)
        real *8 h1s(2,ninterv),h1derss(2,11,ninterv),h1(2)
        real *8 x1,h
c
c 
c   input:
c   ninterv - number of intervals used in equispaced
c       subdivision of current interval
c
c-----------------------------------
c
c       find the subinterval where the point z lives
c
        i=(r-x1)*h +1
        if (i.lt.0) then 
            i = 1
        else if (i.gt.ninterv) then 
            i = ninterv
        endif

        t = (r-x1) - (2*i-1)/h/2
c
c        evaluate the functions h0 and h1 at the point z
c
c
        h0(:)=(((((((((h0derss(:,10,i)*t+
     2      h0derss(:,9,i))*t+h0derss(:,8,i) ) 
     1    * t+h0derss(:,7,i))*t+h0derss(:,6,i))*t+h0derss(:,5,i))*t 
     2    +h0derss(:,4,i))*t+h0derss(:,3,i))*t+h0derss(:,2,i))
     3    *t+h0derss(:,1,i) ) * t + h0s(:,i)
c
        h1(:)=(((((((((h1derss(:,10,i)*t+
     2      h1derss(:,9,i))*t+h1derss(:,8,i) ) 
     1    * t+h1derss(:,7,i))*t+h1derss(:,6,i))*t+h1derss(:,5,i))*t 
     2    +h1derss(:,4,i))*t+h1derss(:,3,i))*t+h1derss(:,2,i))
     3    *t+h1derss(:,1,i) ) * t + h1s(:,i)
c
c
        return
        end
c
c
c
c
c
        subroutine findinte(x,ab,nn,intnum)
        implicit real *8 (a-h,o-z)
        integer intold,ithresh
        dimension ab(2,nn)
c
        data ithresh/10/
c
c       check if the point is on the subinterval as the preceding one
c
        if(intnum .le. 0) goto 2000
        if(intnum .gt. nn) goto 2000
c
        if( (x .ge. ab(1,intnum) ) .and. (x .le. ab(2,intnum) ) ) return
c
 2000 continue
       if(x .lt. ab(1,1)) then
           intnum = 777
           return
       else if(x .gt. ab(2,nn)) then
           intnum = 777
           return
       endif
c
c      the point is not on the same subinterval as the preceding one.
c      if nn is less than ithresh, use direct scan to find the proper 
c      interval
c
       if(nn .gt. ithresh) goto 3000
c
c
        do 2200 j=1,nn
c
           intnum=j
c
        if(ab(2,j) .ge. x) goto 2400
 2200 continue
c
 2400 continue
c
        return
c
 3000 continue
c
c      The point is not on the same subinterval as the preceding one,
c      and nn is greater than ithresh; use bisection to find the proper 
c      interval
c
       i1=1
       i2=nn
       i3=(i1+i2)/2
c
cccc       nsteps=0
       do 3400 i=1,100
c
       if(x .ge. ab(1,i3)) i1=i3
       if(x .le. ab(2,i3)) i2=i3
c
       if(i2 .eq. i1) goto 3600
c
       i3=(i1+i2)/2
 3400 continue
c
 3600 continue
       if(x .lt. ab(1,i3)) i3=i3-1
       if(x .gt. ab(2,i3)) i3=i3+1

       intnum=i3
c       
        return
        end
c
