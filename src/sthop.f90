subroutine sthop

  use mod_global
  implicit none
  
  ! local
  integer, parameter :: nch = 4
  complex (dp) ::  zt, z0, z1
  complex(dp) :: hlp2(nch,nch), hlp1(nch,nch)
  real(dp) :: wc(nch)
  integer :: nf, i_1, ist, i1, i2, i3, i4, m, n, i, j

  ! checkboard breakup for hopping term
  zt = dcmplx(-rt,0.d0)
  do nf = 1,2
     do i_1 = 1,lq/4
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        
        hlp2 = dcmplx(0.d0,0.d0)
        
        n = i1	
        hlp2(1,2) = zt
        hlp2(2,1) = zt
        hlp2(1,4) = zt
        hlp2(4,1) = zt

        n = i2
        hlp2(2,3) = zt
        hlp2(3,2) = zt


        n = i4
        hlp2(4,3) = zt
        hlp2(3,4) = zt

        ! add the chemical potential term
        do i = 1, nch
            hlp2(i,i) = hlp2(i,i) - dcmplx(0.5d0*mu,0.d0) ! 0.5 because every site count 2 times
        end do

        call HermDiag3(nch,nch,hlp2,wc,hlp1)

        do i = 1,nch
           do j = 1,nch
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch
                 z0 = z0 +  hlp1(i,m) * dcmplx(dexp(-dtau *wc(m)),0.d0) * dconjg(hlp1(j,m))
                 z1 = z1 +  hlp1(i,m) * dcmplx(dexp( dtau *wc(m)),0.d0) * dconjg(hlp1(j,m))
              enddo
              urt  (ist,i,j) = z0
              urtm1(ist,i,j) = z1
           enddo
        enddo
     enddo
  enddo

  urt_dn = urt
  urtm1_dn = urtm1

end subroutine sthop
