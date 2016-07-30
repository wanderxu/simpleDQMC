subroutine mmthl(a_up, a_dn)

! perform A * exp(-dtau*T)

  use mod_global
  use matrix_tmp
  implicit none
  
  !arguments:
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_up
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_dn

  ! local
  integer :: nf, i, j, i1, i2, i3, i4, ist

  if (rt.gt.zero) then
      do nf = 2,1,-1
          do i = 1,lq/4
             ist =  i + (nf - 1)*lq/4
             i1 = lthf(i,nf)
             i2 = nnlist(i1,1)
             i3 = nnlist(i1,5)
             i4 = nnlist(i1,2)
             do j = 1,lq
                v1(j) = a_up(j,i1)*urt(ist,1,1) + a_up(j,i2)*urt(ist,2,1) + &
                     &  a_up(j,i3)*urt(ist,3,1) + a_up(j,i4)*urt(ist,4,1)
                v2(j) = a_up(j,i1)*urt(ist,1,2) + a_up(j,i2)*urt(ist,2,2) + &
                     &  a_up(j,i3)*urt(ist,3,2) + a_up(j,i4)*urt(ist,4,2)
                v3(j) = a_up(j,i1)*urt(ist,1,3) + a_up(j,i2)*urt(ist,2,3) + &
                     &  a_up(j,i3)*urt(ist,3,3) + a_up(j,i4)*urt(ist,4,3)
                v4(j) = a_up(j,i1)*urt(ist,1,4) + a_up(j,i2)*urt(ist,2,4) + &
                     &  a_up(j,i3)*urt(ist,3,4) + a_up(j,i4)*urt(ist,4,4)
             enddo
             do j = 1,lq
                a_up(j,i1) = v1(j)
                a_up(j,i2) = v2(j)
                a_up(j,i3) = v3(j)
                a_up(j,i4) = v4(j)
             enddo
          enddo
      enddo
  endif

  if (rt.gt.zero) then
      do nf = 2,1,-1
          do i = 1,lq/4
             ist =  i + (nf - 1)*lq/4
             i1 = lthf(i,nf)
             i2 = nnlist(i1,1)
             i3 = nnlist(i1,5)
             i4 = nnlist(i1,2)
             do j = 1,lq
                v1(j) = a_dn(j,i1)*urt_dn(ist,1,1) + a_dn(j,i2)*urt_dn(ist,2,1) + &
                     &  a_dn(j,i3)*urt_dn(ist,3,1) + a_dn(j,i4)*urt_dn(ist,4,1)
                v2(j) = a_dn(j,i1)*urt_dn(ist,1,2) + a_dn(j,i2)*urt_dn(ist,2,2) + &
                     &  a_dn(j,i3)*urt_dn(ist,3,2) + a_dn(j,i4)*urt_dn(ist,4,2)
                v3(j) = a_dn(j,i1)*urt_dn(ist,1,3) + a_dn(j,i2)*urt_dn(ist,2,3) + &
                     &  a_dn(j,i3)*urt_dn(ist,3,3) + a_dn(j,i4)*urt_dn(ist,4,3)
                v4(j) = a_dn(j,i1)*urt_dn(ist,1,4) + a_dn(j,i2)*urt_dn(ist,2,4) + &
                     &  a_dn(j,i3)*urt_dn(ist,3,4) + a_dn(j,i4)*urt_dn(ist,4,4)
             enddo
             do j = 1,lq
                a_dn(j,i1) = v1(j)
                a_dn(j,i2) = v2(j)
                a_dn(j,i3) = v3(j)
                a_dn(j,i4) = v4(j)
             enddo
          enddo
      enddo
  endif

end subroutine mmthl
