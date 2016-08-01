module ftdqmc_core
  use mod_global
  use matrix_tmp
  use obser
  implicit none

  real(dp), allocatable, dimension(:,:,:), save :: Ust_up, Vst_up, Ust_up_tmp, Vst_up_tmp
  real(dp), allocatable, dimension(:,:,:), save :: Ust_dn, Vst_dn, Ust_dn_tmp, Vst_dn_tmp
  real(dp), allocatable, dimension(:,:), save :: Dst_up, Dst_up_tmp
  real(dp), allocatable, dimension(:,:), save :: Dst_dn, Dst_dn_tmp
  real(dp), dimension(:,:), allocatable, save :: UR_up, VR_up, VL_up, UL_up, Bdtau1_up, grup_tmp
  real(dp), dimension(:,:), allocatable, save :: UR_dn, VR_dn, VL_dn, UL_dn, Bdtau1_dn, grdn_tmp
  real(dp), dimension(:), allocatable, save :: DRvec_up, DLvec_up
  real(dp), dimension(:), allocatable, save :: DRvec_dn, DLvec_dn
  real(dp), dimension(:,:), allocatable, save :: Bt2t1_up, gt0up, g0tup, g00up
  real(dp), dimension(:,:), allocatable, save :: Bt2t1_dn, gt0dn, g0tdn, g00dn


  contains

    subroutine allocate_core
      implicit none
      allocate( Ust_up(ndim,ndim,0:nst) )     ! 1
      allocate( Dst_up(ndim,0:nst) )          ! 2
      allocate( Vst_up(ndim,ndim,0:nst) )     ! 3
      allocate( UR_up(ndim,ndim) )             ! 4
      allocate( DRvec_up(ndim) )               ! 5
      allocate( VR_up(ndim,ndim) )             ! 6
      allocate( VL_up(ndim,ndim) )             ! 7
      allocate( DLvec_up(ndim) )               ! 8
      allocate( UL_up(ndim,ndim) )             ! 9
      allocate( Bdtau1_up(ndim,ndim) )       ! 16
      if(ltau) then
          allocate( Bt2t1_up(ndim,ndim) )       ! 16
          allocate( gt0up(ndim,ndim) )
          allocate( g0tup(ndim,ndim) )
          allocate( g00up(ndim,ndim) )
      end if

      allocate( Ust_up_tmp(ndim,ndim,0:nst) ) ! 17
      allocate( Dst_up_tmp(ndim,0:nst) )      ! 18
      allocate( Vst_up_tmp(ndim,ndim,0:nst) ) ! 19
      allocate( grup_tmp(ndim,ndim) )         ! 20

      allocate( Ust_dn(ndim,ndim,0:nst) )     ! 1
      allocate( Dst_dn(ndim,0:nst) )          ! 2
      allocate( Vst_dn(ndim,ndim,0:nst) )     ! 3
      allocate( UR_dn(ndim,ndim) )             ! 4
      allocate( DRvec_dn(ndim) )               ! 5
      allocate( VR_dn(ndim,ndim) )             ! 6
      allocate( VL_dn(ndim,ndim) )             ! 7
      allocate( DLvec_dn(ndim) )               ! 8
      allocate( UL_dn(ndim,ndim) )             ! 9
      allocate( Bdtau1_dn(ndim,ndim) )       ! 16
      if(ltau) then
          allocate( Bt2t1_dn(ndim,ndim) )       ! 16
          allocate( gt0dn(ndim,ndim) )
          allocate( g0tdn(ndim,ndim) )
          allocate( g00dn(ndim,ndim) )
      end if

      allocate( Ust_dn_tmp(ndim,ndim,0:nst) ) ! 17
      allocate( Dst_dn_tmp(ndim,0:nst) )      ! 18
      allocate( Vst_dn_tmp(ndim,ndim,0:nst) ) ! 19
      allocate( grdn_tmp(ndim,ndim) )         ! 20

    end subroutine allocate_core

    subroutine deallocate_core
      implicit none
      deallocate( grdn_tmp )         ! 20
      deallocate( Vst_dn_tmp )       ! 19
      deallocate( Dst_dn_tmp )       ! 18
      deallocate( Ust_dn_tmp )       ! 17
      if(ltau) then
          deallocate( g00dn )
          deallocate( g0tdn )
          deallocate( gt0dn )
          deallocate( Bt2t1_dn )      ! 16
      end if
      deallocate( Bdtau1_dn )      ! 16
      deallocate( UL_dn )             ! 9
      deallocate( DLvec_dn )          ! 8
      deallocate( VL_dn )             ! 7
      deallocate( VR_dn )             ! 6
      deallocate( DRvec_dn )          ! 5
      deallocate( UR_dn )             ! 4
      deallocate( Vst_dn )         ! 3
      deallocate( Dst_dn )         ! 2
      deallocate( Ust_dn )         ! 1
      deallocate( grup_tmp )         ! 20
      deallocate( Vst_up_tmp )       ! 19
      deallocate( Dst_up_tmp )       ! 18
      deallocate( Ust_up_tmp )       ! 17
      if(ltau) then
          deallocate( g00up )
          deallocate( g0tup )
          deallocate( gt0up )
          deallocate( Bt2t1_up )      ! 16
      end if
      deallocate( Bdtau1_up )      ! 16
      deallocate( UL_up )             ! 9
      deallocate( DLvec_up )          ! 8
      deallocate( VL_up )             ! 7
      deallocate( VR_up )             ! 6
      deallocate( DRvec_up )          ! 5
      deallocate( UR_up )             ! 4
      deallocate( Vst_up )         ! 3
      deallocate( Dst_up )         ! 2
      deallocate( Ust_up )         ! 1
    end subroutine deallocate_core
  
    subroutine ftdqmc_stablize_0b_svd(n)
      ! B( n*tau1, 0 ) = B( n*tau1, (n-1)*tau1 ) * B( (n-1)*tau1, 0 )
      implicit none
      integer, intent(in) :: n

      ! local
      real(dp), allocatable, dimension(:,:) :: Umat1, Umat2, Vmat1, Vmat2
      real(dp), allocatable, dimension(:) :: Dvec1, Dvec2

      allocate( Umat1(ndim,ndim), Umat2(ndim,ndim), Vmat1(ndim,ndim), Vmat2(ndim,ndim) )
      allocate( Dvec1(ndim), Dvec2(ndim) )

      call Bmat_tau( n*nwrap, (n-1)*nwrap+1, Bdtau1_up, Bdtau1_dn )
      
      Umat1(:,:) = Ust_up(:,:,n-1)
      Dvec1(:)   = Dst_up(:,n-1)
      Vmat1(:,:) = Vst_up(:,:,n-1)

      ! Btmp = ( Bdtau1_up * Umat1 ) * Dmat1
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Bdtau1_up,ndim,Umat1,ndim,0.d0,Atmp,ndim)  ! Atmp = Bdtau1_up * Umat1
      call s_d_x_diag_d(ndim,Atmp,Dvec1,Btmp) ! Btmp = Atmp * Dmat1

      call UDVdcmps_R(ndim, ndim, Btmp, Umat2, Dvec2, Vtmp)
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Vtmp,ndim,Vmat1,ndim,0.d0,Vmat2,ndim)  ! Vmat2 = Vtmp * Vmat1
      Ust_up(:,:,n) = Umat2(:,:)
      Dst_up(:,n)   = Dvec2(:)
      Vst_up(:,:,n) = Vmat2(:,:)

      Umat1(:,:) = Ust_dn(:,:,n-1)
      Dvec1(:)   = Dst_dn(:,n-1)
      Vmat1(:,:) = Vst_dn(:,:,n-1)
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Bdtau1_dn,ndim,Umat1,ndim,0.d0,Atmp,ndim)  ! Atmp = Bdtau1_dn * Umat1
      call s_d_x_diag_d(ndim,Atmp,Dvec1,Btmp) ! Btmp = Atmp * Dmat1
      call UDVdcmps_R(ndim, ndim, Btmp, Umat2, Dvec2, Vtmp)
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Vtmp,ndim,Vmat1,ndim,0.d0,Vmat2,ndim)  ! Vmat2 = Vtmp * Vmat1
      Ust_dn(:,:,n) = Umat2(:,:)
      Dst_dn(:,n)   = Dvec2(:)
      Vst_dn(:,:,n) = Vmat2(:,:)

      deallocate( Dvec2, Dvec1 )
      deallocate( Vmat2, Vmat1, Umat2, Umat1 )

    end subroutine ftdqmc_stablize_0b_svd
  
    subroutine ftdqmc_stablize_b0_svd(n)
      ! B( beta, (n-1)*tau1 ) = B( beta, n*tau1 ) * B( n*tau1, (n-1)*tau1 )
      implicit none
      integer, intent(in) :: n

      ! local
      real(dp), allocatable, dimension(:,:) :: Umat1, Umat2, Vmat1, Vmat2
      real(dp), allocatable, dimension(:) :: Dvec1, Dvec2

      allocate( Umat1(ndim,ndim), Umat2(ndim,ndim), Vmat1(ndim,ndim), Vmat2(ndim,ndim) )
      allocate( Dvec1(ndim), Dvec2(ndim) )

      call Bmat_tau( n*nwrap, (n-1)*nwrap+1, Bdtau1_up, Bdtau1_dn )

      Vmat1(:,:) = Vst_up(:,:,n)
      Dvec1(:)   = Dst_up(:,n)
      Umat1(:,:) = Ust_up(:,:,n)

      ! Btmp = Dmat1 * Umat1 * Bdtau1_up
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Umat1,ndim,Bdtau1_up,ndim,0.d0,Atmp,ndim)  ! Atmp = Umat1 * Bdtau1_up
      call s_diag_d_x_d(ndim,Dvec1,Atmp,Btmp) ! Btmp = Dmat1 * Atmp
      call UDVdcmps_R(ndim, ndim, Btmp, Vtmp, Dvec2, Umat2)  ! Btmp = Vtmp * Dmat2 * Umat2
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Vmat1,ndim,Vtmp,ndim,0.d0,Vmat2,ndim)  ! Vmat2 = Vmat1 * Vtmp 
      Vst_up(:,:,n-1) = Vmat2(:,:)
      Dst_up(:,n-1) = Dvec2(:)
      Ust_up(:,:,n-1) = Umat2(:,:)

      Vmat1(:,:) = Vst_dn(:,:,n)
      Dvec1(:)   = Dst_dn(:,n)
      Umat1(:,:) = Ust_dn(:,:,n)
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Umat1,ndim,Bdtau1_dn,ndim,0.d0,Atmp,ndim)  ! Atmp = Umat1 * Bdtau1_dn
      call s_diag_d_x_d(ndim,Dvec1,Atmp,Btmp) ! Btmp = Dmat1 * Atmp
      call UDVdcmps_R(ndim, ndim, Btmp, Vtmp, Dvec2, Umat2)  ! Btmp = Vtmp * Dmat2 * Umat2
      call dgemm('n','n',ndim,ndim,ndim,1.d0,Vmat1,ndim,Vtmp,ndim,0.d0,Vmat2,ndim)  ! Vmat2 = Vmat1 * Vtmp 
      Vst_dn(:,:,n-1) = Vmat2(:,:)
      Dst_dn(:,n-1) = Dvec2(:)
      Ust_dn(:,:,n-1) = Umat2(:,:)

      deallocate( Dvec2, Dvec1 )
      deallocate( Vmat2, Vmat1, Umat2, Umat1 )

    end subroutine ftdqmc_stablize_b0_svd
  
    subroutine ftdqmc_sweep_start
      implicit none
      integer :: n

      ! at tau = 0
      grup(:,:) = Imat(:,:)
      Ust_up(:,:,0) = Imat(:,:)
      Dst_up(:,0)   = Ivec(:)
      Vst_up(:,:,0) = Imat(:,:)

      grdn(:,:) = Imat(:,:)
      Ust_dn(:,:,0) = Imat(:,:)
      Dst_dn(:,0)   = Ivec(:)
      Vst_dn(:,:,0) = Imat(:,:)

      do n = 1, nst
          ! at tau = n * tau1
          call ftdqmc_stablize_0b_svd(n)
      end do
  
      ! at tau = beta
      UR_up(:,:) = Ust_up(:,:,nst)
      DRvec_up(:)= Dst_up(:,nst)
      VR_up(:,:) = Vst_up(:,:,nst)
      call green_equaltime( nst, ndim, UR_up, DRvec_up, VR_up, Imat, Ivec, Imat, grup )

      UR_dn(:,:) = Ust_dn(:,:,nst)
      DRvec_dn(:)= Dst_dn(:,nst)
      VR_dn(:,:) = Vst_dn(:,:,nst)
      call green_equaltime( nst, ndim, UR_dn, DRvec_dn, VR_dn, Imat, Ivec, Imat, grdn )

    end subroutine ftdqmc_sweep_start
  
    subroutine ftdqmc_sweep(lmeasure)
    
      implicit none

      logical, intent(in) :: lmeasure

      ! local variables
      integer :: nt, n, nt1, nt2

      integer, external :: NRanf

      real(dp), allocatable, dimension(:,:) :: Umat1, Vmat1
      real(dp), allocatable, dimension(:) :: Dvec1

      allocate( Umat1(ndim,ndim), Vmat1(ndim,ndim) )
      allocate( Dvec1(ndim) )
  
      ! at tau = beta
      Vst_up(:,:,nst) = Imat(:,:)
      Dst_up(:,nst)   = Ivec(:)
      Ust_up(:,:,nst) = Imat(:,:)

      Vst_dn(:,:,nst) = Imat(:,:)
      Dst_dn(:,nst)   = Ivec(:)
      Ust_dn(:,:,nst) = Imat(:,:)
  
      do nt = ltrot, 1, -1
          if ( mod(nt, nwrap) .eq. 0 .and. (nt/nwrap) .lt. nst ) then
              n = nt/nwrap
              ! at tau = n * tau1
              UR_up(:,:) = Ust_up(:,:,n)
              DRvec_up(:)= Dst_up(:,n)
              VR_up(:,:) = Vst_up(:,:,n)

              UR_dn(:,:) = Ust_dn(:,:,n)
              DRvec_dn(:)= Dst_dn(:,n)
              VR_dn(:,:) = Vst_dn(:,:,n)

              call ftdqmc_stablize_b0_svd(n+1)

              ! for spin up
              UL_up(:,:)  = Ust_up(:,:,n)
              DLvec_up(:) = Dst_up(:,n)  
              VL_up(:,:)  = Vst_up(:,:,n)
              call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp )
              call s_compare_max_d( ndim, grtmp, grup, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
              grup(:,:) = grtmp(:,:)

              ! for spin down
              UL_dn(:,:)  = Ust_dn(:,:,n)
              DLvec_dn(:) = Dst_dn(:,n)  
              VL_dn(:,:)  = Vst_dn(:,:,n)
              call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp )
              call s_compare_max_d( ndim, grtmp, grdn, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
              grdn(:,:) = grtmp(:,:)
          end if
  
          ! obser
          if( lmeasure ) then
             call obser_equaltime(nt)
          end if
  
          !! update
          ! updateu
          if( lupdateu ) then
              call upgradeu( nt, grup, grdn )
              call mmuul  ( grup, grdn, nt )
              call mmuurm1( grup, grdn, nt )
          end if
  
  
          ! wrap H0
          call mmthl  (grup, grdn)
          call mmthrm1(grup, grdn)
  
      end do
  
      ! at tau = 0
      n = 0
      Ust_up(:,:,n) = Imat(:,:)
      Dst_up(:,n)   = Ivec(:)
      Vst_up(:,:,n) = Imat(:,:)

      Ust_dn(:,:,n) = Imat(:,:)
      Dst_dn(:,n)   = Ivec(:)
      Vst_dn(:,:,n) = Imat(:,:)


      if( ltau ) then
          g00up = grup
          gt0up = grup
          g0tup = grup-Imat

          g00dn = grdn
          gt0dn = grdn
          g0tdn = grdn-Imat
      end if
  
      do nt = 1, ltrot, 1

          ! wrap H0
          call mmthr  (grup, grdn)
          call mmthlm1(grup, grdn)

          ! update
          ! updateu
          if( lupdateu ) then
              call mmuur  ( grup, grdn, nt )
              call mmuulm1( grup, grdn, nt )
              if (.not.ltau .or. .not.lmeasure) then
                  call upgradeu( nt, grup, grdn )
              end if
          end if
  
          ! obser
          if( lmeasure .and. .not. ltau ) then
             call obser_equaltime(nt)
          end if
  
          if ( mod(nt, nwrap) .eq. 0 ) then
              n = nt/nwrap
              ! at tau = n * tau1
              VL_up(:,:) = Vst_up(:,:,n)
              DLvec_up(:)= Dst_up(:,n)
              UL_up(:,:) = Ust_up(:,:,n)
              
              VL_dn(:,:) = Vst_dn(:,:,n)
              DLvec_dn(:)= Dst_dn(:,n)
              UL_dn(:,:) = Ust_dn(:,:,n)

              call ftdqmc_stablize_0b_svd(n)

              ! for spin up
              UR_up(:,:)  = Ust_up(:,:,n)
              DRvec_up(:) = Dst_up(:,n)
              VR_up(:,:)  = Vst_up(:,:,n)
              if(.not. ltau ) then
                  call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp )
              else
                  call green_tau(n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, g00up, gt0up, g0tup,  grtmp )
              end if
              call s_compare_max_d( ndim, grtmp, grup, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
              grup(:,:) = grtmp(:,:)

              ! for spin down
              UR_dn(:,:)  = Ust_dn(:,:,n)
              DRvec_dn(:) = Dst_dn(:,n)
              VR_dn(:,:)  = Vst_dn(:,:,n)
              if( .not. ltau ) then
                  call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp )
              else
                  call green_tau(n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, g00dn, gt0dn,  g0tdn,  grtmp )
              end if
              call s_compare_max_d( ndim, grtmp, grdn, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
              grdn(:,:) = grtmp(:,:)

          end if

          ! dyn
          ! g00up g00dn
          ! gt0up gt0dn
          if( ltau .and. lmeasure  ) then

              n = nt/nwrap
              if( mod(nt,nwrap) .eq. 0 ) then

              else
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  call Bmat_tau( nt1, nt2, Bdtau1_up, Bdtau1_dn )

                  ! G(t',0) = B(t',t) * G(t,0)
                  Btmp = gt0up
                  call dgemm('n','n',ndim,ndim,ndim,1.d0,Bdtau1_up,ndim,Btmp,ndim,0.d0,gt0up,ndim)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call InvMatrx_D(ndim,Bdtau1_up)
                  Btmp = g0tup
                  call dgemm('n','n',ndim,ndim,ndim,1.d0,Btmp,ndim,Bdtau1_up,ndim,0.d0,g0tup,ndim)

                  ! G(t',0) = B(t',t) * G(t,0)
                  Btmp = gt0dn
                  call dgemm('n','n',ndim,ndim,ndim,1.d0,Bdtau1_dn,ndim,Btmp,ndim,0.d0,gt0dn,ndim)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call InvMatrx_D(ndim,Bdtau1_dn)
                  Btmp = g0tdn
                  call dgemm('n','n',ndim,ndim,ndim,1.d0,Btmp,ndim,Bdtau1_dn,ndim,0.d0,g0tdn,ndim)

              end if

              call obsert(nt,gt0up,gt0dn,g0tup,g0tdn,grup,grdn)


          end if ! if(ltau) then
  
      end do

      deallocate( Dvec1 )
      deallocate( Vmat1, Umat1 )

    end subroutine ftdqmc_sweep
  
    subroutine green_equaltime( nt, ndm, ure, dre, vre, vle, dle, ule, gtt )
      implicit none
      integer, intent(in) :: nt, ndm
      real(dp), dimension(ndm,ndm), intent(in) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      real(dp), dimension(ndm,ndm), intent(out) :: gtt

      ! local
      integer :: i, j
      real(dp), allocatable, dimension(:,:) :: uutmpinv, uinv_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin

      allocate( uutmpinv(ndm,ndm), uinv_tmp(ndm,ndm) )
      allocate( drmax(ndm), drmin(ndm), dlmax(ndm), dlmin(ndm) )

      ! breakup dre = drmax * drmin
      !         dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)

      ! uutmp = ule*ure
      call dgemm('n','n',ndm,ndm,ndm,1.d0,ule,ndm,ure,ndm,0.d0,uutmp,ndm)  ! uutmp = ule*ure
      ! vvtmp = vre*vle
      call dgemm('n','n',ndm,ndm,ndm,1.d0,vre,ndm,vle,ndm,0.d0,vvtmp,ndm)  ! vvtmp = vre*vle

      uutmpinv = uutmp
      call InvMatrx_D(ndm,uutmpinv)

      !! >> g(t,t)
      ! drmax^-1 * ( ule * ure )^-1 dlmax^-1
      do j = 1, ndm
          do i = 1, ndm
              Atmp(i,j) = uutmpinv(i,j) / ( drmax(i)*dlmax(j) )
          end do
      end do
      ! drmin * ( vre * vle ) * dlmin
      do j = 1, ndm
          do i = 1, ndm
              Btmp(i,j) = vvtmp(i,j) * drmin(i) * dlmin(j)
          end do
      end do

      dvvdtmp(:,:) = Atmp(:,:) + Btmp(:,:)
      call InvMatrx_D(ndm,dvvdtmp)

      uinv_tmp=ule
      call InvMatrx_D(ndm,uinv_tmp)

      call s_v_invd_u( ndm, uinv_tmp, dlmax, dvvdtmp, Btmp )

      uinv_tmp=ure
      call InvMatrx_D(ndm,uinv_tmp)
      
      call s_v_invd_u( ndm, Btmp, drmax, uinv_tmp, gtt )

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( uinv_tmp )
      deallocate( uutmpinv )

    end subroutine green_equaltime

    subroutine green_tau(nt, ndm, ure, dre, vre, vle, dle, ule, g00, gt0, g0t, gtt )
      implicit none
      integer, intent(in) :: nt, ndm
      real(dp), dimension(ndm,ndm), intent(inout) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      real(dp), dimension(ndm,ndm), intent(out) :: g00, gt0, g0t, gtt

      ! local
      integer :: i, j
      real(dp), allocatable, dimension(:,:) :: uutmpinv, uinv_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin

      allocate( uutmpinv(ndm,ndm), uinv_tmp(ndm,ndm) )
      allocate( drmax(ndm), drmin(ndm), dlmax(ndm), dlmin(ndm) )

      ! breakup dre = drmax * drmin
      !         dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)


      ! uutmp = ule*ure
      call dgemm('n','n',ndm,ndm,ndm,1.d0,ule,ndm,ure,ndm,0.d0,uutmp,ndm)  ! uutmp = ule*ure
      ! vvtmp = vre*vle
      call dgemm('n','n',ndm,ndm,ndm,1.d0,vre,ndm,vle,ndm,0.d0,vvtmp,ndm)  ! vvtmp = vre*vle

      uutmpinv = uutmp
      call InvMatrx_D(ndm,uutmpinv)

      !! >> g(t,t)
      ! drmax^-1 * ( ule * ure )^-1 dlmax^-1
      do j = 1, ndm
          do i = 1, ndm
              Atmp(i,j) = uutmpinv(i,j) / ( drmax(i)*dlmax(j) )
          end do
      end do
      ! drmin * ( vre * vle ) * dlmin
      do j = 1, ndm
          do i = 1, ndm
              Btmp(i,j) = vvtmp(i,j) * drmin(i) * dlmin(j)
          end do
      end do

      dvvdtmp(:,:) = Atmp(:,:) + Btmp(:,:)
      call InvMatrx_D(ndm,dvvdtmp)

      uinv_tmp=ule
      call InvMatrx_D(ndm,uinv_tmp)

      call s_v_invd_u( ndm, uinv_tmp, dlmax, dvvdtmp, Btmp )

      uinv_tmp=ure
      call InvMatrx_D(ndm,uinv_tmp)
      
      call s_v_invd_u( ndm, Btmp, drmax, uinv_tmp, gtt )

      !! >> g(t,0)
      call s_v_d_u( ndm, Btmp, drmin, vre, gt0 )

      !! >> g(0,0)
      call InvMatrx_D( ndm, vvtmp )

      ! dlmax^-1 * ( vre * vle )^-1 drmax^-1
      do j = 1, ndm
          do i = 1, ndm
              Atmp(i,j) = vvtmp(i,j) / ( dlmax(i)*drmax(j) )
          end do
      end do
      ! dlmin * ( ule * ure ) * drmin
      do j = 1, ndm
          do i = 1, ndm
              Btmp(i,j) = uutmp(i,j) * dlmin(i) * drmin(j)
          end do
      end do

      dvvdtmp(:,:) = Atmp(:,:) + Btmp(:,:)
      call InvMatrx_D(ndm,dvvdtmp)

      call InvMatrx_D(ndm,vre)

      call s_v_invd_u( ndm, vre, drmax, dvvdtmp, Btmp )

      call InvMatrx_D(ndm,vle)
      call s_v_invd_u( ndm, Btmp, dlmax, vle, g00 )

      !! >> g(0,t)
      call s_v_d_u( ndm, Btmp, dlmin, ule, g0t )
      g0t = -g0t

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( uinv_tmp )
      deallocate( uutmpinv )

    end subroutine green_tau

    subroutine Bmat_tau( nt1, nt2, bmat_up, bmat_dn )
      ! B(tau1,tau2)
      implicit none
      integer, intent(in) :: nt1, nt2
      real(dp), dimension(ndim,ndim), intent(out) :: bmat_up
      real(dp), dimension(ndim,ndim), intent(out) :: bmat_dn

      ! local
      integer :: nt

      bmat_up(:,:) = Imat(:,:)
      bmat_dn(:,:) = Imat(:,:)
      do nt = nt2, nt1
          call mmthr(bmat_up,bmat_dn)
          if( lupdateu ) then
            call mmuur(bmat_up, bmat_dn, nt )
          end if
      end do
    end subroutine Bmat_tau
end module ftdqmc_core
