subroutine preq
  use mod_global
  use obser
  implicit none

  integer :: i, j
  complex(dp) :: expiqr, spipi

  ! calculate  the S(pi,pi)
  ! S(pi,pi) = 1/L^2 \sum_ij <  1/2 (ni_up -ni_dn) *  1/2 ( nj_up - nj_dn ) >  exp( i * (pi,pi) * r(i-j) )
  ! normalize
  spin_corrlt(:,:) = spin_corrlt(:,:) / dcmplx( dble(nobs), 0.d0 )
  spipi = czero
  do i = 1, lq
      do j = 1, lq
          expiqr = exp( dcmplx( 0.d0, ( pi*dble(list(j,1)-list(i,1)) + pi*dble( list(j,2)-list(i,2) ) ) ) )
          spipi  = spipi + spin_corrlt(j,i) * expiqr
      end do
  end do

  spipi  = spipi  / dcmplx( dble( 4 * lq * lq), 0.d0 )  ! 1/4 comes from 1/2 spin

  open (unit=40,file='spin_corrlt.bin',status='unknown', action="write", position="append")

  write(40, '((2e16.8,2x))') spipi

  close(40)

  ! calculate obs_bin
  obs_bin(:) = obs_bin(:) / dcmplx( dble(nobs), 0.d0 )
  open (unit=90,file='ener1.bin',status='unknown', action="write", position="append")
  write(90, '(3(e16.8,4x))') dble(obs_bin(1))/dble(ndim), dble(obs_bin(2))/dble(ndim), dble(obs_bin(3))/dble(ndim)
  close(90)


end subroutine preq
