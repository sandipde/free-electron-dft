subroutine del2(V,P)
use kind_param
use global_var
implicit none
real(R8B),dimension(ndim_x,ndim_y,ndim_z)::V,P
integer(I4B)::i,j,k
!allocate (V(nx,ny,nz))
!allocate (P(nx,ny,nz))
P=0.0
do i=2,ndim_x-1
 do j=2,ndim_y-1
  do k=2,ndim_z-1
   
   P(i,j,k)=(V(i,j,k+1)+V(i,j,k-1)-2.0*V(i,j,k))/delta_z**2+(V(i,j+1,k)+V(i,j-1,k)-2.0*V(i,j,k))/delta_y**2+(V(i+1,j,k)+V(i-1,j,k)-2.0*V(i,j,k))/delta_x**2
  enddo
 enddo
enddo
P=-P/12.56637061
return
end subroutine del2
