subroutine sli
  use mod_global

  implicit none

  ! local
  integer :: ncount, nx, ny, n, i, j, iq, ix, iy, nk, imj_nx, imj_ny, imj
  real(dp) :: ri(2), qvec(2), rr

  logical :: ltest
  integer :: nf, nc, i0, nc1, nc2

  integer, external :: npbc

  real(dp), dimension(:,:), allocatable :: rlist

  list = 0; invlist = 0; nnlist = 0

  ncount = 0
  do nx = 1,l
  do ny = 1,l
     ncount = ncount + 1
     list(ncount,1) = nx
     list(ncount,2) = ny
     invlist(nx,ny) = ncount
  enddo
  enddo

  do n = 1,lq
     nx = list(n,1)
     ny = list(n,2)
     nnlist(n,0) = invlist( nx , ny )
     nnlist(n,1) = invlist( npbc(nx+1,l) , ny )
     nnlist(n,2) = invlist( nx , npbc(ny+1,l) )
     nnlist(n,3) = invlist( npbc(nx-1,l) , ny )
     nnlist(n,4) = invlist( nx , npbc(ny-1,l) )
     nnlist(n,5) = invlist( npbc(nx+1,l) , npbc(ny+1,l) )
     nnlist(n,6) = invlist( npbc(nx-1,l) , npbc(ny+1,l) )
     nnlist(n,7) = invlist( npbc(nx-1,l) , npbc(ny-1,l) )
     nnlist(n,8) = invlist( npbc(nx+1,l) , npbc(ny-1,l) )
  enddo
  
  ! latt_imj
  do j = 1, lq
      do i = 1, lq
          imj_nx = npbc( list(i,1) - list(j,1), l )
          imj_ny = npbc( list(i,2) - list(j,2), l )
          latt_imj(i,j) = invlist( imj_nx, imj_ny )
      end do
  end do

  ! latt_listk
  nk = 0
  do j = 0, l-1
      do i = 0, l-1
          nk = nk+1
          listk(nk,1) = i-l/2
          listk(nk,2) = j-l/2
      end do
  end do
  if( nk .ne. lq ) then
      stop " Error in set listk "
  end if

  ! set zexpiqr
  do iq = 1, lq
      qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
      do i = 1, lq
          ri = dble(list(i,1))*a1_p + dble(list(i,2))*a2_p
          zexpiqr(i,iq) = exp( dcmplx( 0.d0,  qvec(1)*ri(1) + qvec(2)*ri(2) ) )
      end do
  end do

  
  ! set sites with special coordinate for decompostion of kinetic matrix exp(-dtau*T)
  nc1 = 0
  nc2 = 0
  do i = 1,lq
     ix = list(i,1)
     iy = list(i,2)
     if (mod(ix,2).eq.0 .and. mod(iy,2).eq.0 ) then
        nc1 = nc1 + 1
        lthf(nc1,1) = i
     endif
     if (mod(ix,2).ne.0 .and. mod(iy,2).ne.0 ) then
        nc2 = nc2 + 1
        lthf(nc2,2) = i
     endif
  enddo
  if( l .gt. 1 ) then
      if (nc1.ne.lq/4 .or. nc2  .ne. lq/4 ) then
         write(6,*) 'error 1'
         stop
      endif
  end if

end subroutine sli

integer function npbc(nr,l)
  implicit none
  integer, intent(in) :: nr
  integer, intent(in) :: l
  npbc = nr
  if (nr.gt.l) npbc = nr - l
  if (nr.lt.1) npbc = nr + l
end function npbc
