subroutine double_grid_trans(charge_den)
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
integer(I4B)::x_origin,y_origin,z_origin
real (R8B)::charge_den(ndim_x,ndim_y,ndim_z)
!_____________________________________________________________________
do i=-ndim_x+1,ndim_x
	d_x(i)=i*delta_x
enddo
do i=-ndim_y+1,ndim_y
	d_y(i)=i*delta_y
enddo
do i=-ndim_z+1,ndim_z
	d_z(i)=i*delta_z
enddo
original_x_start=(ndim_x)/2
original_y_start=(ndim_y)/2
original_z_start=(ndim_z)/2
!write(*,*) original_x_start,original_y_start,original_z_start
double_charge_den=0.0
!write(*,*)"111"
do i=-original_x_start+1,original_x_start
  do j=-original_y_start+1,original_y_start
    do k=-original_z_start+1,original_z_start
!     write(*,*)"ijk=",i+original_x_start,j+original_y_start,k+original_z_start
      double_charge_den(i,j,k)=charge_den(i+original_x_start,j+original_y_start,k+original_z_start)
    enddo
  enddo
enddo
!write(*,*)"1234"
!open (2,file='double_charge_den.dat')
!do i=-ndim_x,ndim_x
!  do j=-ndim_y,ndim_y
!    write(2,200)d_x(i),d_y(j),d_z(ndim_z),double_charge_den(i,j,0)
!  enddo 
!    write(2,*)"  "
!enddo
!200 Format (4(F15.7,2X))
return
end subroutine double_grid_trans
