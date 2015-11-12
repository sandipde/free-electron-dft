!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________

subroutine charge_density()
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
kappa  =  (0.944863d0)/ang2bohr
kappa  =  a0/2.0
gauss_const=dsqrt((kappa*kappa/pi)**3.0)
!write(*,*)"g_const=",gauss_const
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
      rho_up(i,j,k)=gauss_const*dexp(-kappa*kappa*(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)))
    enddo
   enddo
enddo
open (1,file='charge_den.dat')
do i=1,ndim_x
  do j=1,ndim_y
    write(1,200)xi(i),yi(j),zi(ndim_z),rho_up(i,j,ndim_z/2)
  enddo 
    write(1,*)"  "
enddo
200 Format (4(F15.7,2X))
end subroutine charge_density

!=================================================================================================

subroutine gen_charged_sphere(rmax)
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
real(R8B)::r,rmax,vol_sphere
!rmax=5.0
vol_sphere=4.0/3.0*pi*rmax*rmax*rmax
 rho_up=0.0
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
	r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
        if (r.le.rmax)then
     !        charge_den(i,j,k)=1.0/(vol_sphere)*4.0/3.0*pi*r**3/(vol_sphere)
        rho_up(i,j,k)=1.0/(vol_sphere)
	endif 
    enddo
   enddo
enddo
! charge_den=0.0
open (1,file='charge_den.dat')
do i=1,ndim_x
  do j=1,ndim_y
   do k=1,ndim_z
       write(1,*) dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)),rho_up(i,j,k)
   enddo
  enddo
  enddo
200 Format (4(F15.7,2X))
end subroutine gen_charged_sphere
