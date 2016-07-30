module mod_global
  integer, parameter :: dp = 8
  integer, parameter :: fout = 50
  integer, parameter :: finp = 51
  real(dp), parameter :: zero = 0.d0
  complex(dp), parameter :: cone = dcmplx( 1.d0, 0.d0 )
  complex(dp), parameter :: czero = dcmplx( 0.d0, 0.d0 )
  complex(dp), parameter :: chalf = dcmplx( 0.5d0, 0.d0 )
  complex(dp), parameter :: cquarter = dcmplx( 0.25d0, 0.d0 )
  real(dp), save :: pi = dacos(-1.d0)
 
  ! lattice
  integer, save :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
  real(dp), save :: b1_p(2), b2_p(2)
  integer, save :: ndim
  integer, save :: l
  integer, save :: lq
  integer, save :: nfam
  integer, save :: lfam
  integer, allocatable, dimension(:,:), save :: list
  integer, allocatable, dimension(:,:), save :: invlist
  integer, allocatable, dimension(:,:), save :: nnlist
  integer, allocatable, dimension(:,:), save :: list_plaq
  integer, allocatable, dimension(:,:), save :: ltpf
  integer, allocatable, dimension(:,:), save :: lthf

  integer, allocatable, dimension(:,:), save :: latt_imj
  integer, allocatable, dimension(:,:), save :: listk

  ! model
  integer, save :: ne
  real(dp), parameter :: rt = 1.d0
  real(dp), save :: beta
  real(dp), save :: mu
  real(dp), save :: rhub

  ! dqmc relative
  logical, save :: lupdateu
  integer, save :: iseed
  real(dp), save :: dtau
  integer, save :: ltrot
  complex(dp), allocatable, dimension(:,:), save :: Imat
  complex(dp), allocatable, dimension(:,:), save :: grup, grdn, grupc, grdnc
  real(dp), allocatable, dimension(:), save :: Ivec

  ! cal. control
  integer, save :: nbin
  integer, save :: nsweep
  integer, save :: nst
  integer, save :: nwrap
  logical, save :: lwarnup
  integer, save :: nwarnup
  integer, save :: n_outconf_pace 

  ! for dynamical
  logical, save :: ltau
  complex(dp), allocatable, dimension(:,:), save :: zexpiqr

  ! for obs
  integer, save :: nobs
  integer, save :: obs_segment_len 
  complex(dp), dimension(10), save :: main_obs, mpi_main_obs


  ! for DQMC
  integer, dimension(-1:1,1), save :: nflipl
  complex(dp), dimension(-1:1), save :: xsigma_u_up, xsigma_u_dn
  complex(dp), dimension(-1:1,1), save :: delta_u_up, delta_u_dn 

  complex(dp), allocatable, dimension(:,:,:), save :: urt, urtm1
  complex(dp), allocatable, dimension(:,:,:), save :: urt_dn, urtm1_dn

  integer, allocatable, dimension(:,:), save :: nsigl_u
  real(dp), save :: max_wrap_error, max_wrap_error_tmp 

  contains

  subroutine make_tables
    implicit none

    integer :: i
    logical :: exists

    ! default parameters
    l    = 4
    beta = 4
    dtau = 0.05d0
    mu   = 0.d0 ! default is half filling
    rhub = 2.0d0
    nwrap = 10
    n_outconf_pace = 1

    ltau = .false.

    nsweep = 200
    nbin = 5
    obs_segment_len = 10

    lupdateu = .false.
    
    ! read parameters
    exists = .false.
    inquire (file = 'ftdqmc.in', exist = exists)
    if ( exists .eqv. .true. ) then
        open( unit=finp, file='ftdqmc.in', status='unknown' )
        read(finp,*) rhub, l, beta, dtau, nwrap, nsweep, nbin, ltau
    else
        write(fout,'(a)') ' No ftdqmc.in found, start simulation with default parameters '
    end if

    ! tune parameters
    if( rhub .gt. -0.001d0 ) lupdateu = .true.
    lq = l*l
    lfam = max(lq/2,1)
    nfam = 1
    ndim = lq   ! the dimension of matrix inside determinant
    ltrot = nint( beta / dtau )
    nst = ltrot / nwrap

   	a1_p(1) = 1 ; a1_p(2) =  0
    a2_p(1) = 0 ; a2_p(2) =  1
    L1_p = l*a1_p
    L2_p = l*a2_p

    b1_p(1) = 2.d0*pi/dble(l) ; b1_p(2) = 0.d0
    b2_p(1) = 0.d0            ; b2_p(2) = 2.d0*pi/dble(l)

    ! allocate tables
    allocate( list(lq,2) )
    allocate( invlist(l,l) )
    allocate( nnlist(lq, 0:8) )
    allocate( list_plaq(lq, 1:5) )
    allocate( ltpf(max(lq/2,1), 4) )
    allocate( lthf(max(lq/4,1), 2) )

    allocate( latt_imj(lq,lq) )
    allocate( listk(lq,2) )

    allocate( zexpiqr(lq,lq) )

    allocate( urt(max(lq/2,1),4,4), urtm1(max(lq/2,1),4,4) )
    allocate( urt_dn(max(lq/2,1),4,4), urtm1_dn(max(lq/2,1),4,4) )

    allocate(grup(ndim,ndim), grdn(ndim,ndim), grupc(ndim,ndim), grdnc(ndim,ndim))

    allocate( Imat(ndim,ndim) )
    allocate( Ivec(ndim) )
    Imat = czero
    do i = 1, ndim
        Imat(i,i) = cone
        Ivec(i) = 1.d0
    end do

  end subroutine make_tables

  subroutine deallocate_tables
    deallocate( Ivec, Imat )
    deallocate( grdnc, grupc, grdn, grup )
    deallocate( urtm1_dn, urt_dn )
    deallocate( urtm1, urt )

    deallocate( zexpiqr )

    deallocate( listk, latt_imj )
    deallocate( lthf, ltpf, list_plaq, nnlist, invlist, list )
  end subroutine deallocate_tables

end module mod_global
