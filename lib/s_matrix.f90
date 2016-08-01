  subroutine s_compare_max_z( ndim, Amat, Bmat, max_diff )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  compare two martrix Amat and Bmat, and return the maximum difference
  !%  Input: ndim, Amat, Bmat
  !%  Output: max_diff
    implicit none

    integer, intent(in) :: ndim
    complex(8), dimension(ndim,ndim), intent(in) :: Amat, Bmat
    real(8), intent(out) :: max_diff
    
    ! local
    integer :: i, j
    real(8) :: max_tmp

    max_diff = 0.d0
    do j = 1, ndim
        do i = 1, ndim
            max_tmp = abs( Amat(i,j) - Bmat(i,j) )
            if( max_tmp .gt. max_diff ) then
                max_diff = max_tmp
            end if
        end do
    end do
  end subroutine s_compare_max_z

  subroutine s_compare_max_d( ndim, Amat, Bmat, max_diff )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  compare two martrix Amat and Bmat, and return the maximum difference
  !%  Input: ndim, Amat, Bmat
  !%  Output: max_diff
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim,ndim), intent(in) :: Amat, Bmat
    real(8), intent(out) :: max_diff
    
    ! local
    integer :: i, j
    real(8) :: max_tmp

    max_diff = 0.d0
    do j = 1, ndim
        do i = 1, ndim
            max_tmp = dabs( Amat(i,j) - Bmat(i,j) )
            if( max_tmp .gt. max_diff ) then
                max_diff = max_tmp
            end if
        end do
    end do
  end subroutine s_compare_max_d

  subroutine s_z_x_diag_d( ndim, Amat, dvec, Bmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform the product of a full matrix and a diagnoal matrix
  !%  Input: ndim, Amat, dvec
  !%  Output: Bmat
    implicit none
    integer, intent(in) :: ndim
    complex(8), dimension(ndim,ndim), intent(in) :: Amat
    real(8), dimension(ndim), intent(in) :: dvec
    complex(8), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dcmplx( dvec(i), 0.d0 ) * Amat(j,i)
        end do
    end do

  end subroutine s_z_x_diag_d

  subroutine s_d_x_diag_d( ndim, Amat, dvec, Bmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform the product of a full matrix and a diagnoal matrix
  !%  Input: ndim, Amat, dvec
  !%  Output: Bmat
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim,ndim), intent(in) :: Amat
    real(8), dimension(ndim), intent(in) :: dvec
    real(8), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dvec(i) * Amat(j,i)
        end do
    end do

  end subroutine s_d_x_diag_d

  subroutine s_diag_d_x_z( ndim, dvec, Amat, Bmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform the product of a full matrix and a diagnoal matrix
  !%  Input: ndim, Amat, dvec
  !%  Output: Bmat
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: dvec
    complex(8), dimension(ndim,ndim), intent(in) :: Amat
    complex(8), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dcmplx( dvec(j), 0.d0 ) * Amat(j,i)
        end do
    end do

  end subroutine s_diag_d_x_z

  subroutine s_diag_d_x_d( ndim, dvec, Amat, Bmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform the product of a full matrix and a diagnoal matrix
  !%  Input: ndim, Amat, dvec
  !%  Output: Bmat
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: dvec
    real(8), dimension(ndim,ndim), intent(in) :: Amat
    real(8), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dvec(j) * Amat(j,i)
        end do
    end do

  end subroutine s_diag_d_x_d

  subroutine s_diag_dvd( ndim, dvecr, Amat, dvecl, Bmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform diagnoal matrix * full matrix * diagnoal matrix
  !%  Input: ndim, dvecr, Amat, dvecl
  !%  Output: Bmat
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: dvecr, dvecl
    complex(8), dimension(ndim,ndim), intent(in) :: Amat
    complex(8), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dcmplx( dvecr(j)*dvecl(i), 0.d0 ) * Amat(j,i)
        end do
    end do

  end subroutine s_diag_dvd

  subroutine s_v_invd_u( ndim, vmat, dvec, umat, zmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform full matrix * diagnoal matrix ^ -1 * full matrix
  !%  Input: ndim, vmat, dvec, umat
  !%  Output: zmat
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: dvec
    real(8), dimension(ndim,ndim), intent(in) :: vmat, umat
    real(8), dimension(ndim,ndim), intent(out) :: zmat
    
    ! local
    real(8) :: ztmp
    integer :: i, j, k
    do i = 1, ndim
        do j = 1, ndim
            ztmp = 0.d0
            do k = 1, ndim
                ztmp = ztmp + vmat(j,k)*umat(k,i) / dvec(k)
            end do
            zmat(j,i) = ztmp
        end do
    end do

  end subroutine s_v_invd_u

  subroutine s_v_d_u( ndim, vmat, dvec, umat, zmat )
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform full matrix * diagnoal matrix * full matrix
  !%  Input: ndim, vmat, dvec, umat
  !%  Output: zmat
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: dvec
    real(8), dimension(ndim,ndim), intent(in) :: vmat, umat
    real(8), dimension(ndim,ndim), intent(out) :: zmat
    
    ! local
    real(8) :: ztmp
    integer :: i, j, k
    do i = 1, ndim
        do j = 1, ndim
            ztmp = 0.d0
            do k = 1, ndim
                ztmp = ztmp + vmat(j,k)*umat(k,i) * dvec(k)
            end do
            zmat(j,i) = ztmp
        end do
    end do

  end subroutine s_v_d_u

  subroutine s_dvec_min_max(ndim,dvec,dmax,dmin)
  !%  Written by Xiao Yan Xu (wanderxu@gmail.com)
  !%  Date: July 30, 2016
  !%  perform the decompostion of a vector, dvec = dmax * dmin, to force all numbers greater than one in dmax, all numbers less than one in dmin
  !%  Input: ndim, dvec
  !%  Output: dmax, dmin
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: dvec
    real(8), dimension(ndim), intent(out) :: dmax, dmin

    ! local
    integer :: i

    do i = 1, ndim
        if( dvec(i) .ge. 1.d0 ) then
            dmin(i) = 1.d0
            dmax(i) = dvec(i)
        end if
        if( dvec(i) .le. 1.d0 ) then
            dmax(i) = 1.d0
            dmin(i) = dvec(i)
        end if
    end do

  end subroutine s_dvec_min_max
