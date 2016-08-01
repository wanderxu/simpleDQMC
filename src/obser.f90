module obser
  use mod_global

  real(dp), allocatable, dimension(:,:), save :: spin_corrlt

  real(dp), save :: obs_bin(10)

  real(dp), allocatable, dimension(:,:), save :: gtau

  contains

  subroutine allocate_obs
    implicit none
    allocate( spin_corrlt(lq,lq) )
    if(ltau) then
        allocate( gtau(ndim,ltrot) )
    end if
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    if(ltau) then
        deallocate( gtau)
    end if
    deallocate( spin_corrlt )
  end subroutine deallocate_obs

  subroutine obser_init
    implicit none
    nobs = 0
    spin_corrlt(:,:) = 0.d0
    obs_bin(:) = 0.d0
    if(ltau) then
        gtau(:,:) = 0.d0
    end if

  end subroutine obser_init

  subroutine obser_equaltime(nt)
    implicit none
    integer,intent(in) :: nt

    ! local 
    integer :: i, j, i1, i2, i3, i4, i_1, nf, ist
    real(dp) :: szsz_tmp, zkint, zedoub, zne

    nobs = nobs + 1

    !grup (i,j) = < c_i c^+_j >
    !grupc (i,j) = < c^+_i c_j >

    ! get grupc
    do i = 1, ndim
        do j = 1, ndim
            grupc(j,i) = - grup(i,j)
        end do
        grupc(i,i) = grupc(i,i) + 1
    end do

    ! get grdn and grdnc
    do i = 1, ndim
        do j = 1, ndim
            grdnc(j,i) = - grdn(i,j)
        end do
        grdnc(i,i) = grdnc(i,i) + 1
    end do

    ! zne
    zne = 0.d0
    do i = 1, ndim
        zne = zne + grupc(i,i) + grdnc(i,i)
    end do
    obs_bin(1) = obs_bin(1) + zne

    ! measure kinetic energy
    zkint = 0.d0
    IF ( l .gt. 1 ) THEN
    do nf = 1,2
       do i_1 = 1,lq/4
          ist = i_1 + (nf - 1)*lq/4
          i1 = lthf(i_1,nf)
          i2 = nnlist(i1,1)
          i3 = nnlist(i1,5)
          i4 = nnlist(i1,2)
          
          zkint = zkint  +  grupc(i1,i2) + &
                            grupc(i2,i1)
          zkint = zkint  +  grupc(i1,i4) + &
                            grupc(i4,i1) 

          zkint = zkint  +  grdnc(i1,i2) + &
                            grdnc(i2,i1)
          zkint = zkint  +  grdnc(i1,i4) + &
                            grdnc(i4,i1) 

          zkint = zkint  +  grupc(i2,i3) + &
                            grupc(i3,i2) 

          zkint = zkint  +  grdnc(i2,i3) + &
                            grdnc(i3,i2)

          zkint = zkint  +  grupc(i4,i3) + &
                            grupc(i3,i4) 

          zkint = zkint  +  grdnc(i4,i3) + &
                            grdnc(i3,i4) 
       end do
    end do
    ELSE
        zkint = zkint - 4.d0*rt * ( grupc(1,1) + grdnc(1,1) )
    ENDIF
    obs_bin(2) = obs_bin(2) + zkint*(-rt)

    zedoub = 0.d0
    do i = 1, lq
        zedoub = zedoub + grupc(i,i)*grdnc(i,i)
    end do

    obs_bin(3) = obs_bin(3) + zedoub

    ! measure fermion spin sz sz correlation
    !   ( ni_up -ni_dn ) * ( nj_up - nj_dn ) 
    ! = ni_up * nj_up + ni_dn * nj_dn - ni_dn * nj_up - ni_up * nj_dn
    ! = grupc(i,i)*grupc(j,j) + grupc(i,j)*grup(i,j) + grdnc(i,i)*grdnc(j,j) + grdnc(i,j)*grdn(i,j) - grdnc(i,i)*grupc(j,j) - grupc(i,i)*grdnc(j,j)
    !                                                  
    do i = 1, lq
        do j = 1, lq
            szsz_tmp = grupc(i,i)*grupc(j,j) + grupc(i,j)*grup(i,j) + &
                       grdnc(i,i)*grdnc(j,j) + grdnc(i,j)*grdn(i,j) - &
                       grdnc(i,i)*grupc(j,j) - grupc(i,i)*grdnc(j,j)
            spin_corrlt(j,i) = spin_corrlt(j,i) + szsz_tmp
        end do
    end do

  end subroutine obser_equaltime

  subroutine obsert(nt, grt0_up, grt0_dn, gr0t_up, gr0t_dn, grtt_up, grtt_dn)
    implicit none
    integer, intent(in) :: nt
    real(dp), dimension(ndim,ndim), intent(in) :: grt0_up, grt0_dn, gr0t_up, gr0t_dn, grtt_up, grtt_dn

    ! local 
    integer :: i, j, imj, iax, imx, jax, jmx

    do j = 1, lq
        jax = nnlist(j,1)
        jmx = nnlist(j,3)
        do i = 1, lq
            imj  = latt_imj(i,j)
            iax = nnlist(i,1) ! i+x
            imx = nnlist(i,3) ! i-x

            gtau(imj,nt) = gtau(imj,nt) + (grt0_up(i,j)+grt0_dn(i,j))*0.5d0
        end do
    end do
  end subroutine obsert

end module obser
