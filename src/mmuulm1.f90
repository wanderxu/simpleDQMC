subroutine mmuulm1(a_up, a_dn, ntau )

! perform A * exp(-V(c))

  use mod_global
  use matrix_tmp
  implicit none
  
  ! arguments:
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_up
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_dn
  integer, intent(in) :: ntau
  
  ! local
  integer :: nl, j
  
  do j  = 1,ndim
     do nl = 1,ndim
        a_up(nl,j) = a_up(nl,j) /  xsigma_u_up(nsigl_u(j,ntau))
        a_dn(nl,j) = a_dn(nl,j) /  xsigma_u_dn(nsigl_u(j,ntau))
     enddo
  enddo
  
end subroutine mmuulm1
