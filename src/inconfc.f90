subroutine inconfc

    use mod_global
    implicit none

    ! local
    integer :: iseed0, i, nt
    real(dp) :: x

    real(dp), external :: Ranf

    allocate( nsigl_u(lq,ltrot) )
        
    open( unit=30,file='confin', status='unknown' )

	read(30,*) iseed0
	if (iseed0.eq.0) then
       ! start from scratch
       lwarnup = .true.
       write(fout,'(a)') ' start from scratch, need warnup '
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
       do nt = 1,ltrot
          do i  = 1,lq
              read(30,*) nsigl_u(i,nt)
          enddo
       enddo
    end if

    close(30)

end subroutine inconfc
