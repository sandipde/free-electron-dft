subroutine write_igenvecup_for_vmd(Eigen_vec)
use kind_param
use global_var
implicit none
integer::i,j,k,state
real(R8B)::volume,Eigen_vec(IWRSZ)
volume =hamilton_dim*delta_x*delta_y*delta_z
100 Format (3(F8.5,2X))
!_____________________________________________________________________________________
! GENERATING THE COMMIN FORMAT FILE FOR VMD
open(8,file='../output_for_vmd/vmdformat')
write(8,*) "QD"
write(8,*) 1
write(8,100) ndim_x*delta_x,0, 0
write(8,100) 0,ndim_y*delta_y,0
write(8,100) 0, 0,ndim_z*delta_z
write(8,*) 1
write(8,*) "Direct"
write(8,100) ndim_x*delta_x*0.5, ndim_y*delta_y*0.5, ndim_z*delta_z*0.5
write(8,*)""
write(8,*) ndim_x,ndim_y,ndim_z
 close(8)
!_______________________________________________________________________________________
!_________GNERATING EIGEN VECTOR FILES.EACH FILE CONTAINS 6 VECTORS COLUMN WISE._________
if (NUMEMAX>=1)then
open(9,file='../output_for_vmd/1st6_eigen_vecup.dat')
do i=1,hamilton_dim
write(9,300) (volume*Eigen_vec(i+j*hamilton_dim),j=0,5)
enddo
endif
if (NUMEMAX>6)then
open(10,file='../output_for_vmd/7to12_eigen_vecup.dat')
do i=1,hamilton_dim
write(10,300) (volume*Eigen_vec(i+j*hamilton_dim),j=6,11)
enddo
endif
if (NUMEMAX>12)then
open(11,file='../output_for_vmd/13to18_eigen_vecup.dat')
do i=1,hamilton_dim
write(11,300) (volume*Eigen_vec(i+j*hamilton_dim),j=12,17)
enddo
endif
if (NUMEMAX>17)then
open(12,file='../output_for_vmd/19to24_eigen_vecup.dat')
do i=1,hamilton_dim
write(12,300) (volume*Eigen_vec(i+j*hamilton_dim),j=18,23)
enddo
endif
if (NUMEMAX>24)then
open(13,file='../output_for_vmd/25to30_eigen_vecup.dat')
do i=1,hamilton_dim
write(13,300) (volume*Eigen_vec(i+j*hamilton_dim),j=24,29)
enddo
endif
300 Format (6(E12.6,3X))
!____________________________________________________________________________________
RETURN
end subroutine write_igenvecup_for_vmd
