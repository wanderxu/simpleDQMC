subroutine upgradeu(ntau, green_up, green_dn)

  use mod_global
  use matrix_tmp

  implicit none

  !arguments
  integer,intent(in) :: ntau
  complex(dp), intent(inout), dimension(ndim,ndim) :: green_up, green_dn

  !local
  complex(dp) ::  ratioup, ratiodn, ratiotot, del44_up, del44_dn
  integer :: i1, nl, nl1, nl2, nrflip
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  real(dp), external :: Ranf

  accm  = 0.d0
  do i1 = 1,lq
     nrflip = 1

     del44_up   =  delta_u_up( nsigl_u(i1,ntau), nrflip )
     ratioup = dcmplx(1.d0,0.d0) + del44_up * ( cone - green_up(i1,i1) )

     del44_dn   =  delta_u_dn( nsigl_u(i1,ntau), nrflip )
     ratiodn = dcmplx(1.d0,0.d0) + del44_dn * ( cone - green_dn(i1,i1) )

     ratiotot = (ratioup*ratiodn)

     ratio_re = dble( ratiotot )

     ratio_re_abs = ratio_re
     if (ratio_re .lt. 0.d0 )  ratio_re_abs = - ratio_re 

     random = Ranf(iseed)
     if ( ratio_re_abs .gt. random ) then

        accm  = accm + 1.d0

        weight = dsqrt(dble(ratiotot*dconjg(ratiotot)))

         
        ! update greep_up
        do nl = 1, ndim
            u1(nl) = green_up(nl,i1)/ratioup
            v1(nl) = green_up(i1,nl)
        end do
        v1(i1) = v1(i1) - cone  ! note the sign
        do nl = 1, ndim
            v1(nl) = del44_up * v1(nl)
        end do

        do nl2 = 1,ndim
        do nl1 = 1,ndim
           green_up(nl1,nl2) = green_up(nl1,nl2) + u1(nl1)*v1(nl2)
        enddo
        enddo

        ! update greep_dn
        do nl = 1, ndim
            u1(nl) = green_dn(nl,i1)/ratiodn
            v1(nl) = green_dn(i1,nl)
        end do
        v1(i1) = v1(i1) - cone  ! note the sign
        do nl = 1, ndim
            v1(nl) = del44_dn * v1(nl)
        end do

        do nl2 = 1,ndim
        do nl1 = 1,ndim
           green_dn(nl1,nl2) = green_dn(nl1,nl2) + u1(nl1)*v1(nl2)
        enddo
        enddo

        ! flip
        nsigl_u(i1,ntau) =  nflipl(nsigl_u(i1,ntau), nrflip)
        
     endif
  enddo
  main_obs(1) = main_obs(1) + dcmplx( accm, dble(lq) )
end subroutine upgradeu
