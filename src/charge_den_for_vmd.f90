subroutine charge_den_for_vmd()
use kind_param
use global_var
implicit none
integer::i,j,k
real(R8B)::volume
volume =hamilton_dim*delta_x*delta_y*delta_z
100 Format (3(F8.5,2X))
!_____________________________________________________________________________________
! GENERATING THE Up_spin charge_den FOR VMD
open(101,file='../chargeden_for_vmd/up_spin.chgcar')
write(101,*) "QD"
write(101,*) 1
write(101,100) ndim_x*delta_x,0, 0
write(101,100) 0,ndim_x*delta_x,0
write(101,100) 0, 0,ndim_z*delta_z
write(101,*) 1
write(101,*) "Direct"
write(101,100) ndim_x*delta_x*0.5, ndim_y*delta_y*0.5, ndim_z*delta_z*0.5
write(101,*)""
write(101,*) ndim_x,ndim_y,ndim_z
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
   write(101,*)volume*rho_up(i,j,k)
  enddo
 enddo
enddo
 close(101)
!_______________________________________________________________________________________

!_____________________________________________________________________________________
! GENERATING THE down_spin charge_den FOR VMD
open(101,file='../chargeden_for_vmd/dn_spin.chgcar')
write(101,*) "QD"
write(101,*) 1
write(101,100) ndim_x*delta_x,0, 0
write(101,100) 0,ndim_x*delta_x,0
write(101,100) 0, 0,ndim_z*delta_z
write(101,*) 1
write(101,*) "Direct"
write(101,100) ndim_x*delta_x*0.5, ndim_y*delta_y*0.5, ndim_z*delta_z*0.5
write(101,*)""
write(101,*) ndim_x,ndim_y,ndim_z
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
   write(101,*)volume*rho_dn(i,j,k)
  enddo
 enddo
enddo
 close(101)
!_______________________________________________________________________________________

!_____________________________________________________________________________________
! GENERATING THE Up_spin+dn_spin charge_den FOR VMD
open(101,file='../chargeden_for_vmd/tot_charge_den.chgcar')
write(101,*) "QD"
write(101,*) 1
write(101,100) ndim_x*delta_x,0, 0
write(101,100) 0,ndim_x*delta_x,0
write(101,100) 0, 0,ndim_z*delta_z
write(101,*) 1
write(101,*) "Direct"
write(101,100) ndim_x*delta_x*0.5, ndim_y*delta_y*0.5, ndim_z*delta_z*0.5
write(101,*)""
write(101,*) ndim_x,ndim_y,ndim_z
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
   write(101,*)volume*(rho_up(i,j,k)+rho_dn(i,j,k))
  enddo
 enddo
enddo
 close(101)
!_______________________________________________________________________________________

!_____________________________________________________________________________________
! GENERATING THE Up_spin-dn_spin charge_den FOR VMD
open(101,file='../chargeden_for_vmd/spin_den.chgcar')
write(101,*) "QD"
write(101,*) 1
write(101,100) ndim_x*delta_x,0, 0
write(101,100) 0,ndim_x*delta_x,0
write(101,100) 0, 0,ndim_z*delta_z
write(101,*) 1
write(101,*) "Direct"
write(101,100) ndim_x*delta_x*0.5, ndim_y*delta_y*0.5, ndim_z*delta_z*0.5
write(101,*)""
write(101,*) ndim_x,ndim_y,ndim_z
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
   write(101,*)volume*(rho_up(i,j,k)-rho_dn(i,j,k))
  enddo
 enddo
enddo
 close(101)
!_______________________________________________________________________________________

return
end subroutine charge_den_for_vmd
