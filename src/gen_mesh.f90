!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________
subroutine gen_mesh()
use kind_param
use global_var
implicit none
integer(I4B)::i,x_origin,y_origin,z_origin

x_bound=(ndim_x)*delta_x/2.0
y_bound=(ndim_y)*delta_y/2.0
z_bound=(ndim_z)*delta_z/2.0
Lx=2*x_bound
Ly=2*y_bound
Lz=2*z_bound
 cell_vol=Lx*Ly*Lz
write(50,*)"****************************MESH PARAMETERS*****************************"
write(50,100)ndim_x,ndim_y,ndim_z
write(50,'(10X,a35,2x,f5.3)')"delta_x=",delta_x
write(50,'(10X,a35,2x,f5.3)')"delta_y=",delta_y
write(50,'(10X,a35,2x,f5.3)')"delta_z=",delta_z
write(50,*)"     Lx =",Lx,"Bohr=",Lx*bohr2ang,"Angstrom"
write(50,*)"     Ly =",Ly ,"Bohr=",Ly*bohr2ang,"Angstrom"
write(50,*)"     Lz =",Lz ,"Bohr=",Lz*bohr2ang,"Angstrom"
100 Format (20X,'Mesh dimension::',3(I4,3X))
do i=1,ndim_x
	xi(i)=-x_bound+(i-1)*delta_x
!	if (xi(i)==0.0) x_origin=i
enddo
do i=1,ndim_y
	yi(i)=-y_bound+(i-1)*delta_y
!	if (yi(i)==0.0) y_origin=i
enddo
do i=1,ndim_z
	zi(i)=-z_bound+(i-1)*delta_z
!	if (zi(i)==0.0)  z_origin=i
enddo
!write(*,*) "ORIGIN OF THE MESH IS AT::(",x_origin,y_origin,z_origin,")"
write(*,*)"______________________________________________________________________________"
return
end subroutine gen_mesh
