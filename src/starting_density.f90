!=======================================================================================
!SUBROUTINE FOR CREATING INITIAL CHARGE DENSITY
!DEPENDING UPON THE INPUT THIS WILL CREATE OR READ CHARGE DENSITIES
!PROGRAM BY : SANDIP DE  DATE 17.02.09
!========================================================================================

subroutine starting_density()
use kind_param
use global_var
implicit none
integer(i4B)::spin,occupied_state,i,j,k
real(R8B),allocatable::rho(:,:,:)!,vh(:,:,:),igenE(:)
real(R8B)::intg,norm

allocate(rho(ndim_x,ndim_y,ndim_z))

if (ans=="y")then

!--------------------CONSTRUCTING RHO FIRST TIME----------------------------------
 write(50,*)"==================== CONSTRUCTING RHO FIRST TIME ========================="
 call solve_3dhamilton(ext_pot,igenE_up,Eigen_vec_up,0)
 call modpsi2(Eigen_vec_up,norm)
! write(*,*)"norm =",norm
 igenE_dn=igenE_up  ! for 1st time only
 call fillup_electrons(occupied_stateup,occupied_statedn)
do spin=1,2 			! spin=1 :: up 
 if (spin==1)then
   occupied_state=occupied_stateup
!   rho=rho_up
 else
   occupied_state=occupied_statedn
 !  rho=rho_dn
 endif
 call gen_rho(spin,occupied_state,Eigen_vec_up,rho)
 if (spin==1) then
   rho_up=rho
 else
   rho_dn=rho
 endif
 call int_rho(rho,intg)
enddo
 call check_rho(0)
write(50,*)"====================== CONSTRUCTION SUCCESSFUL================================"
deallocate(rho)

else
!-------------------------------------------------------------------------------------------------
if (cont=="y") then
write(50,*)"====================CONTINUING UNFINISHED JOB======================="
 open (unit=125,file='../back_up/emergency_rho_up.dat',form='unformatted',status='old')
 open (unit=126,file='../back_up/emergency_rho_dn.dat',form='unformatted',status='old')
else
 if (start_from_backup=="y") then
     if (last_iter==0)then   
 	write(50,*)"=====================READING starting CHARGE DENSITIES FROM FILE======================="
 	open (unit=125,file='../back_up/startrho_up.dat',form='unformatted',status='old')
 	open (unit=126,file='../back_up/startrho_dn.dat',form='unformatted',status='old')
     else
       write(50,*)"=====================READING CHARGE DENSITIES FROM FILE======================="
       open (unit=125,file='../back_up/rho_up.dat',form='unformatted',status='old')
       open (unit=126,file='../back_up/rho_dn.dat',form='unformatted',status='old')
     endif
   endif
endif
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
   read(125)rho_up(i,j,k)
   read(126)rho_dn(i,j,k)
  enddo
 enddo
enddo
!open (60,file='../info/energy.info',status='old')
!write(*,*)"last backup iter?"
!read(*,*)last_iter
!last_iter=last_iter+2
!do i=1,last_iter
!read(60,*)
!enddo
endif

return
end subroutine starting_density
