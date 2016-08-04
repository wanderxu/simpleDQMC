module mod_global
! this module contains:
!  1. shared data
!     constants, tables for lattice, parameters for model, DQMC related matrices, control parameters
!  2. subroutines
!     make_tables : read input, allocate tables
!     deallocate_tables : deallocate tables allocated in make_tables

  ! constants
  integer, parameter :: dp = 8
  integer, parameter :: fout = 50
  integer, parameter :: finp = 51
  real(dp), save :: pi = dacos(-1.d0)
 
  ! lattice
  integer, save :: a1_p(2), a2_p(2)  ! primitive cell vector
  integer, save :: L1_p(2), L2_p(2)  ! lattice size vector
  real(dp), save :: b1_p(2), b2_p(2) ! least interval of k points

  integer, save :: ndim  ! dimension of B matrix, here equals l^2
  integer, save :: l     ! length of lattice
  integer, save :: lq    ! l^2
  integer, save :: nfam  ! number of family for Ising field
  integer, allocatable, dimension(:,:), save :: list      ! coordinate of sites
  integer, allocatable, dimension(:,:), save :: invlist   ! index of sites
  integer, allocatable, dimension(:,:), save :: nnlist    ! neighbors of sites
  integer, allocatable, dimension(:,:), save :: lthf      ! sites with special coordinate for decompostion of kinetic matrix exp(-dtau*T)

  integer, allocatable, dimension(:,:), save :: latt_imj  ! index of ri-rj
  integer, allocatable, dimension(:,:), save :: listk     ! coordinate of k points

  ! model
  real(dp), parameter :: rt = 1.d0   ! hopping amplitude
  real(dp), save :: beta             ! inverse of temperature
  real(dp), save :: mu               ! chemical potential
  real(dp), save :: rhub             ! Hubbard U

  ! dqmc relative
  logical, save :: lupdateu          ! control update
  integer, save :: iseed             ! seed for random number generator
  real(dp), save :: dtau             ! trotter decompostion step
  integer, save :: ltrot             ! total number of time slice after trotter decompostion
  real(dp), allocatable, dimension(:,:), save :: Imat   ! Identity matrix
  real(dp), allocatable, dimension(:,:), save :: grup, grdn, grupc, grdnc   ! green functions
  real(dp), allocatable, dimension(:), save :: Ivec        ! vector will all elements one

  ! cal. control
  integer, save :: nbin         ! number of bins
  integer, save :: nsweep       ! number of sweep
  integer, save :: nwrap        ! numeric stablization pace
  integer, save :: nst          ! number of numberic stablizations
  logical, save :: lwarnup      ! control warmup
  integer, save :: nwarnup      ! number of sweep for warmup
  integer, save :: n_outconf_pace  ! the pace of configuration output

  ! for dynamical
  logical, save :: ltau       ! control dynamic measurement
  complex(dp), allocatable, dimension(:,:), save :: zexpiqr  ! table for exp(i*q*r)

  ! for obs
  integer, save :: nobs             ! number of obser in one bin
  real(dp), dimension(10), save :: main_obs !  acculators for some observation in main program, like accept rate


  ! for DQMC
  integer, dimension(-1:1), save :: nflipl      ! for flipping auxilary field
  real(dp), dimension(-1:1), save :: xsigma_u_up, xsigma_u_dn  ! interaction part matrix element e^(-V(c)), depend on filed
  real(dp), dimension(-1:1), save :: delta_u_up, delta_u_dn  ! ratio of interaction part matrix element for different fileds

  real(dp), allocatable, dimension(:,:,:), save :: urt, urtm1       ! part of e^(-dtau*T) and e^(dtau*T) for up spin flavor
  real(dp), allocatable, dimension(:,:,:), save :: urt_dn, urtm1_dn ! part of e^(-dtau*T) and e^(dtau*T) for dn spin flavor

  integer, allocatable, dimension(:,:), save :: nsigl_u                ! auxilary field
  real(dp), save :: max_wrap_error, max_wrap_error_tmp                 ! error of progating green function

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

    lupdateu = .false.
    
    ! read parameters
    exists = .false.
    inquire (file = 'ftdqmc.in', exist = exists)
    if ( exists .eqv. .true. ) then
        open( unit=finp, file='ftdqmc.in', status='unknown' )
        read(finp,*) rhub, l, beta, dtau, nwrap, nsweep, nbin, ltau
        close(finp)
    else
        write(fout,'(a)') ' No ftdqmc.in found, start simulation with default parameters '
    end if

    ! tune parameters
    if( rhub .gt. 0.d0 ) lupdateu = .true.
    lq = l*l
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
    allocate( lthf(max(lq/4,1), 2) )

    allocate( latt_imj(lq,lq) )
    allocate( listk(lq,2) )

    allocate( zexpiqr(lq,lq) )

    allocate( urt(max(lq/2,1),4,4), urtm1(max(lq/2,1),4,4) )
    allocate( urt_dn(max(lq/2,1),4,4), urtm1_dn(max(lq/2,1),4,4) )

    allocate(grup(ndim,ndim), grdn(ndim,ndim), grupc(ndim,ndim), grdnc(ndim,ndim))

    allocate( Imat(ndim,ndim) )
    allocate( Ivec(ndim) )
    Imat = 0.d0
    do i = 1, ndim
        Imat(i,i) = 1.d0
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
    deallocate( lthf, nnlist, invlist, list )
  end subroutine deallocate_tables

end module mod_global
