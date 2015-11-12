
!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________
subroutine read_input()
use kind_param
use global_var
implicit none
open(1,file="../input/mesh.in")
read(1,*);read(1,*) ndim_x,ndim_y,ndim_z
read(1,*);read(1,*) delta_x,delta_y,delta_z
 close (1)
open(77,file='../input/main.in')
read(77,*)n_eup					!up e no
read(77,*)n_edn					!dn e no
read(77,*)deg_limit				!max diff in energy for degenarate levels
read(77,*)mix_para				!charge density mixing parameter
read(77,*)conv_limit				!energy conv_limit
read(77,*)max_iter				!max iter no
!read(77,*)max_coeff                             ! number of expansion coeff is to be consider in finite difference expression of second derivative
max_coeff=6
read(77,*)backup_iter                           ! no of iteration after which backup is to be created
read(77,*)max_change_allowed                    !Max change allowed in energy in to consecutive iterations
read(77,*)ans					!if ans=y initial charge density will be created else read from file
read(77,*)fstop                                 !if fstop=y stop after creating initial density
read(77,*)start_from_backup			!start from last sceduled backup?
read(77,*)cont					!y to continue from a last saved point
read(77,*)last_iter				!last iteration after which backup was created

if (ans=='n' .and.fstop=='y') then
 write(50,'(10X,"ERROR in main input !!!")')
 stop
endif 
!if (start_from_backup==cont) then
! write(50,'(10X,"ERROR in main input (not sure from where to start)!!!")')
! stop
!endif 
write(50,'(10X,"YOU HAVE",I4,2X,"UP-SPIN & ",I4,2X,"DOWN-SPIN ELECTRONS")') n_eup,n_edn
write(50,'("BACKUP WILL BE CREATED AFTER EACH ",I4," ITERATIONS")')backup_iter
write(50,'("MAXIMUM ENERGY CHANGE ALLOWED IN ONE ITERATION ",F10.5)')max_change_allowed
!write(50,*)"CONTINUING FROM A LAST SAVED POINT"
return

end subroutine read_input
