subroutine salph

  use mod_global
  implicit none

  ! local 
  real(dp) :: alpha_u

  write(fout,*)
  write(fout,'(a)') ' In salph: '

  nflipl(-1, 1) =  1
  nflipl( 1, 1) = -1
  
  ! onsite
  ! spin up, positive coupling
  alpha_u = acosh(dexp(0.5d0*dtau*rhub))
  xsigma_u_up(-1) = dcmplx( dexp(  alpha_u * ( -1.d0 ) ), 0.d0 )  ! -1.d0 is the field, 0.5 is from fermion spin s=1/2 
  xsigma_u_up( 1) = dcmplx( dexp(  alpha_u * (  1.d0 ) ), 0.d0 )  !  1.d0 is the field, 0.5 is from fermion spin s=1/2
  ! spin down, negative coupling
  xsigma_u_dn(-1) = dcmplx( dexp( -alpha_u * ( -1.d0 ) ), 0.d0 )  ! -1.d0 is the field, 0.5 is from fermion spin s=1/2
  xsigma_u_dn( 1) = dcmplx( dexp( -alpha_u * (  1.d0 ) ), 0.d0 )  !  1.d0 is the field, 0.5 is from fermion spin s=1/2

  delta_u_up(-1,1) = ( xsigma_u_up( 1) / xsigma_u_up(-1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from -1 to  1
  delta_u_up( 1,1) = ( xsigma_u_up(-1) / xsigma_u_up( 1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from  1 to -1

  delta_u_dn(-1,1) = ( xsigma_u_dn( 1) / xsigma_u_dn(-1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from -1 to  1
  delta_u_dn( 1,1) = ( xsigma_u_dn(-1) / xsigma_u_dn( 1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from  1 to -1

  write(fout,'(a,e18.8)') ' alpha_u = ', alpha_u

end subroutine salph
