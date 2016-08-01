subroutine prtau
  use mod_global
  use obser
  implicit none

  character(40) :: filek

  interface
     subroutine fourier_trans_tau(gr,filek)
       real(kind=8), dimension(:,:) :: gr
       character (40) :: filek
     end subroutine fourier_trans_tau
  end interface

  gtau= gtau/dble(nsweep)

  filek = "gtau.bin"
  call fourier_trans_tau(gtau,filek)
end subroutine prtau

subroutine fourier_trans_tau(gr,filek)
  use mod_global
  implicit none
  real(dp), dimension(:,:) :: gr
  integer :: imj, nt, nk
  character (40) :: filek
  real(dp) :: xk_p(2)
  complex(dp), allocatable , dimension(:,:) :: gk

  allocate (gk(lq,ltrot))

  gk = dcmplx(0.d0,0.d0)
  do imj = 1,lq
     do nt = 1,ltrot
        do nk = 1,lq
           gk(nk,nt) = gk(nk,nt) +  gr(imj,nt)/zexpiqr(imj,nk)
        enddo
     enddo
  enddo
  gk = gk/dcmplx(dble(lq),0.d0)

  open (unit=20,file=filek,status='unknown', action="write", position="append")
  do nk = 1,lq
     xk_p = dble(listk(nk,1))*b1_p + dble(listk(nk,2))*b2_p
     write(20,*) xk_p(1), xk_p(2)
     do nt = 1,ltrot
         write(20,*) gk(nk,nt)
     enddo
  enddo
  close(20)

  deallocate (gk)
end subroutine fourier_trans_tau
