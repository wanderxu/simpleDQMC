subroutine outconfc

    use mod_global
    implicit none

    ! local
    integer :: i, nt

    open( unit=35,file='confout', status='unknown' )

    write(35,*) iseed

    do nt = 1,ltrot
       do i  = 1,lq
           write(35,'(i2)') nsigl_u(i,nt)
       enddo
    enddo

    close(35)

end subroutine outconfc
