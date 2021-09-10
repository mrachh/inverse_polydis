        subroutine hank106_wrap(r,h0s,h1s,w,lw,intnum,n)
        implicit real *8 (a-h,o-z)
        dimension r(n)
        complex *16 h0s(n),h1s(n)

        ifexpon = 1

        do i=1,n
        call hank106_r_h01(r(i),h0s(i),h1s(i),ifexpon,intnum,w,lw)
        enddo

        return
        end
