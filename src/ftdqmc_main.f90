program ftdqmc_main

!>  finite temperature determinant quantum Monte Carlo package
!>  Here is an example for solving Hubbard model on square lattice
!> 
!>  Hamiltionian:
!>     H = -t \sum_<ij>s c_is^dagger c_js + h.c. + U \sum_i ( n_iup - 1/2) * ( n_idn - 1/2 )
!>   
!>  Partition function
!>     Z = \sum_c w(c),   c is configuration
!>
!>     with w(c) = det( 1 + Bup(M)*Bup(M-1)*....*Bup(2)Bup(1) ) * det( 1 + Bdn(M)*Bdn(M-1)*....*Bdn(2)Bdn(1) )
!> 
!>     Bs(tau) = exp(s*alpha_u*Diag(S_tau)) * exp( -dtau * T )
!>     

  use mod_global
  use matrix_tmp
  use ftdqmc_core
  implicit none

  ! local
  integer :: nbc, nsw
  real(dp) :: start_time, end_time, time1, time2

  call cpu_time(start_time)

  open( unit=fout, file='ftdqmc.out', status='unknown' )

  main_obs(:) = 0.d0

  call ftdqmc_initial

  ! read in parameters and allocate tables
  call make_tables

  ! set lists for lattice
  call sli

  call ftdqmc_initial_print
 
  ! prepare for the DQMC
  call salph   ! set phase of Ising field
  call inconfc ! initial configuration
  call sthop   ! compute exp(-dtau*T)
  


  call allocate_matrix_tmp ! allocate tmp data for matrix operations
  call allocate_core       ! allocate data for sweep
  call allocate_obs        ! allocate data for obserables

  max_wrap_error = 0.d0

  call ftdqmc_sweep_start

  write(fout,'(a)') ' ftdqmc_sweep_start done '
  write(fout,*)

  ! warnup
  if( lwarnup ) then
      ! set nwarnup
      nwarnup = nint(4*lq*beta)  ! usually, 4*lq*beta sweeps of warnup is enough.
      if(rhub.le.0.d0) nwarnup = 0
      write(fout,'(a,i8)') ' nwarnup = ', nwarnup
      do nsw = 1, nwarnup
          call ftdqmc_sweep(lmeasure=.false.)
      end do
      write(fout,'(a)') ' warmup done '
      write(fout,*)
  end if

  call cpu_time(time1)
  do nbc =  1, nbin

      call obser_init

      do nsw = 1, nsweep

          call ftdqmc_sweep(lmeasure=.true.)

      end do

      call preq            ! output equaltime measurement data to bins
      if(ltau) call prtau  ! output dynamical measurement data to bins

      ! output configuration control
      if( nbc .eq. 1 )  then
          call cpu_time(time2)
          n_outconf_pace = nint( dble( 3600 * 12 ) / ( time2-time1 ) ) ! set output configuration pace, default is every 12 hours
          if( n_outconf_pace .lt. 1 ) n_outconf_pace = 1
          write(fout,'(a,e16.8,a)') ' time for 1 bin: ', time2-time1, ' s'
          write(fout,'(a,i12)') ' n_out_conf_pace = ', n_outconf_pace
          write(fout,*)
      end if

      if( n_outconf_pace .lt. nbin/3 ) then  ! output configuraton every n_outconf_pace bins, if it is less than nbin/3
          if( mod(nbc,n_outconf_pace) .eq. 0 ) then
              call outconfc
          end if
      else if( mod( nbc, max(nbin/3,1) ) .eq. 0 ) then ! output configuration every nbin/3 bins if n_outconf_pace is too large
          call outconfc
      end if

      write( fout, '(i5,a,i5,a)' ) nbc, '  /', nbin, '   finished '

  end do

  write(fout,*)
  write(fout, '(a,e16.8)') ' max_wrap_error = ', max_wrap_error

  call outconfc

  write(fout,*)
  if(lupdateu)  write(fout,'(a,e16.8)') ' >>> accep_u  = ', main_obs(1)/main_obs(2)


  call deallocate_core
  call deallocate_matrix_tmp

  call deallocate_tables

  call cpu_time(end_time)
  write(fout,*)
  write(fout,'(a,f10.2,a)') ' >>> Total time spent:', end_time-start_time, 's'
  write(fout,*)
  write(fout,'(a)') ' The simulation done !!! '
  write(fout,*)
  write(fout,'(a)') '        o         o    '
  write(fout,'(a)') '       o o       o o   '
  write(fout,'(a)') '       o o       o o   '
  write(fout,'(a)') '        o         o    '
  write(fout,'(a)') '       o o       o o   '
  write(fout,'(a)') '       o o       o o   '
  write(fout,'(a)') '        o         o    '

  close(fout)

end program ftdqmc_main
