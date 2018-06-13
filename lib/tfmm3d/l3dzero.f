c
c
        subroutine l3dzero(mpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       ... set multipole to zero
c
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
        return
        end
c
c
c
