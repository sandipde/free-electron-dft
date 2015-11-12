subroutine analytic_pot(rmax,exact_pot)
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
real(R8B)::exact_pot(ndim_x,ndim_y,ndim_z),rmax,r,vol_sphere
vol_sphere=4.0*pi*rmax**3/3.0
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
     r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
     if (r <= rmax) then
        exact_pot(i,j,k)=(3.0*rmax*rmax-r*r)/(8.0*pi*rmax**3)
      else
        exact_pot(i,j,k)=1.0/(4.0*pi*r)
     endif
    enddo
  enddo
enddo
open (3,file='analytic_pot.dat')
do i=1,ndim_x
  do j=1,ndim_y
   do k=1,ndim_z
       write(3,*) dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)),exact_pot(i,j,k)
   enddo
  enddo
  enddo
return
end subroutine analytic_pot
