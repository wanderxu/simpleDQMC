subroutine ftdqmc_initial
  use mod_global

  integer :: system_time

  !=========================================================
  !%% inital the seed for pseudo random number generator   $
  !---------------------------------------------------------
  call system_clock(system_time)
  iseed = abs( system_time - 1014748364 )*(pi/dsqrt(19.d0))

  ! print head
  write(fout,'(a)') ' ===================================================================================='
  write(fout,*)
  write(fout,'(a)') '        The finite temperature determinant quantum monte carlo (DQMC) package '
  write(fout,*)
  write(fout,'(a)') '            FFFF   TTTTT   DDD      QQQ     M   M      CCCC                    '
  write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
  write(fout,'(a)') '            FFFF     T     D   D   Q   Q   M M M M    C                        '
  write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
  write(fout,'(a)') '            F        T     DDD      QQQ   M   M   M    CCCC                    '
  write(fout,'(a)') '                                      \                                       '
  write(fout,*)
  write(fout,*)
  write(fout,'(a)') ' written by Xiao Yan Xu ( wanderxu@gmail.com )                                '
  write(fout,*)
  write(fout,'(a)') ' history: '
  write(fout,*)
  write(fout,'(a)') '     30/07/2016,  version 1.0  '
  write(fout,*)
  write(fout,'(a)') ' ------------------------------------------------------------------------------------'
  write(fout,*)
  write(fout,'(a)') ' >>> The simulation start running now ! '
end subroutine ftdqmc_initial

subroutine ftdqmc_initial_print
  use mod_global
  implicit none

  integer :: i, j

  write(fout,*)
  write(fout,'(a)')' --------------------------------- '
  write(fout,'(a)')' We will simulate with parameters  '
  write(fout,'(a)')' --------------------------------- '
  write(fout,*)
  write(fout,'(a,f6.2)')    ' t      = ', rt
  write(fout,'(a,f6.2)')    ' U      = ', rhub
  write(fout,'(a,f7.3)')    ' mu     = ', mu
  write(fout,'(a,f6.2)')    ' beta   = ', beta
  write(fout,'(a,i4)')      ' L      = ', l
  write(fout,'(a,f7.3)')    ' dtau   = ', dtau
  write(fout,'(a,i6)')      ' ltrot  = ', ltrot
  write(fout,'(a,i6)')      ' nwrap  = ', nwrap
  write(fout,'(a,i6)')      ' nsweep = ', nsweep
  write(fout,'(a,i6)')      ' nbin   = ', nbin
  write(fout,*)  'lupdateu = ', lupdateu
  write(fout,*)  'ltau = ', ltau

  write(fout,*)
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)')' The lattice sites list '
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)') '   i     list(i,:) '
  do i = 1, lq
      write(fout,'(i6,2i4)') i,  list(i,:)
  end do

end subroutine ftdqmc_initial_print
