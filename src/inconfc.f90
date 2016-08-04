subroutine inconfc

    use mod_global
    implicit none

    ! local
    integer :: iseed0, i, nt
    real(dp) :: x
    logical :: exists

    real(dp), external :: Ranf

    allocate( nsigl_u(lq,ltrot) )
        
    exists = .false.
    inquire (file = 'confin', exist = exists)
    if ( exists .eqv. .true. ) then
        open( unit=30,file='confin', status='unknown' )
	    read(30,*) iseed0
    else
        iseed0 = 0
        write(fout,'(a)') ' No confin file found '
    end if
	if (iseed0.eq.0) then
       if( exists .eqv. .true. ) close(30)
       ! start from scratch
       lwarnup = .true.
       write(fout,'(a)') ' start from scratch, need warnup '
       write(fout,*)
	   do nt = 1,ltrot
           do i  = 1,lq
	           x = Ranf(iseed)
	           nsigl_u(i,nt) = 1
	           if (x.gt.0.5) nsigl_u(i,nt) = -1
           enddo
       enddo
	else
       iseed = iseed0
       ! read configurations from confin
       lwarnup = .false.
       write(fout,'(a)') ' start from old configuration, do not need warnup '
       write(fout,*)
       do nt = 1,ltrot
          do i  = 1,lq
              read(30,*) nsigl_u(i,nt)
          enddo
       enddo
       close(30)
    end if

end subroutine inconfc
