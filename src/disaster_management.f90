subroutine disaster_management(iter)
use kind_param
use global_var
implicit none
integer(I4B)::iter,i,j,k
real(R8B)::time1,time2
if (((iter-(iter/backup_iter)*backup_iter)==0) .OR. iter==999)then     !write after each backup_iter iteration
open (80,file= '../info/disaster_management.info')
write(80,*)"DISASTER MANAGEMENT STARTED !!!!!"
 call cpu_time(time1)
!open (35,file='../output/Eigen_vec_up.dat')
!open (45,file='../output/Eigen_vec_dn.dat')
!do i=1,IWRSZ
! write(35,'(E20.10)')Eigen_vec_up
! write(45,'(E20.10)')Eigen_vec_dn
!enddo
 call write_igenvecup_for_vmd(Eigen_vec_up)
 call write_igenvecdn_for_vmd(Eigen_vec_dn)
if (iter==999)then
open (unit=125,file='../back_up/emergency_rho_up.dat',form='unformatted')
open (unit=126,file='../back_up/emergency_rho_dn.dat',form='unformatted')
else
if (iter==0)then
open (unit=125,file='../back_up/startrho_up.dat',form='unformatted')
open (unit=126,file='../back_up/startrho_dn.dat',form='unformatted')
else
open (unit=125,file='../back_up/rho_up.dat',form='unformatted')
open (unit=126,file='../back_up/rho_dn.dat',form='unformatted')
endif
endif
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
   write(125)rho_up(i,j,k)
   write(126)rho_dn(i,j,k)
  enddo
 enddo
enddo

if (iter==999) then
  write (80,*) "Doing an Emergency backup..."
endif
write(80,*)"BACKUP CREATED AFTER ITERATION NO >",iter
 call cpu_time(time2)
write(80,*)"time taken =",(time2-time1),"secs"
write(*,*)
endif
return
end subroutine disaster_management
