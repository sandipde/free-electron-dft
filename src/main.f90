! SANDIP DE
! Date:01.12.08
!______________________________________________________________________________
program main
use kind_param					! here the parameters are defined
use global_var					! All global variables are defined here
implicit none
integer(i4b)::i,j,k,iter,start_time(3),end_time(3),spin,occupied_state,em_backup         
real(R8B)::norm,intg,tot_E,temp_E,start_iter,end_iter,sec_start,sec_end
!_______________________________________________________________________________



!			GENERAL INFORMATIONS
!====================================================================================
 open (50,file='../info/run.info')                      ! all details are written in this file
 write(50,*)"========================================================================="
 call itime(start_time)
 call cpu_time(sec_start)
 write(50,'(10X,"PROGRAM STARTS AT",3I10)')start_time
 !======================================================================================
 write(*,'("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
 write(*,*)
 write(*,*)"TO SEE WHAT IS GOING ON DO THE FOLLOWING >>>>>>"
 write(*,*)"cd ../info/"
 write(*,*)"cat <filename>.info (where filename = run,energy,davidson,fermi_occupancy,eigen,rho_check)"
 write(*,*)"BACKUPS ARE CREATED IN ../back_up/ FOLDER"
 write(*,*)"FINAL WAVE FUNCTIONS ARE IN THE ../output_for_vmd FOLDER.."
 write(*,*)
 write(*,'("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
!=========================================================================================================






!-------------------------INITIALIZATION--------------------------------------------
 call read_input()                        ! Read all the inputs
 n_e=n_eup+n_edn			  ! totall_electron=up_el+dn_el
 call dvd_setup()			  ! Davidson setup 
 kT=0.01				  ! This does not have any effect still now
 call nonzero_sparse_counter()            ! this calculates the requried dim of the compressed hamilton
 include 'allocate'                       ! all allocations are done here
 call gen_mesh()			  ! generate Mesh
 call ext_potential()			  ! generate External potential based on inputs
 call calc_rs()                           ! Calculate Rs of the system
!---------------------------------------------------------------------------------

last_iter=0
if (ans=="n" )then
 if(cont=='y')then
  open (60,file='../info/energy.info',status='old')
  write(*,*)"last iter?"
  read(*,*)last_iter
  do i=1,last_iter+2
  read(60,*)
  enddo
 endif
  if(start_from_backup=='y')then
  open (60,file='../info/energy.info',status='old')
!  write(*,*)"last backup iter?"
!  read(*,*)last_iter
  do i=1,last_iter+2
  read(60,*)
  enddo
 endif
!write(60,'(2X,"----------------------------------------------------------------------------------------------------------------------------")')
else 
open (60,file='../info/energy.info')
write(60,'(4X,"ITER NO",10X,"KINETIC ENERGY",10X,"E_ext",10X,"E_hart",13X,"E_xc",13X,"TOTAL ENERGY",6X,"TIME (SECS)")')
write(60,'(2X,"=============================================================================================================================")')
endif


!=======================================STARTING CALCULATIONS=================================
 call starting_density()                  ! The initial density is constructed here
 
 call disaster_management(0)              ! A backup in case any thing goes wrong

 if(fstop=="y") stop
 call gen_kinetic_hamilton()      	  ! calculate del2 which will be operated on wave fuction to find kinetic energy
!==============================================================================================



!--------------------------------------------------------------------------------------------------------
write(50,*)"SELF CONSISTENCY STARTING >>>>>>>>"
write(50,'(10X,"CONVERGENCE ENEGRY DIFFERENCE LIMIT = ",E10.4)')conv_limit
write(50,'(10X,"MAXIMUM NO OF ITERATION = ",I6)')max_iter
write(50,'("^^^^^^^^^^^^MIXING ",F5.2," % OF NEW RHO WITH OLD ONE ^^^^^^^^^^^^^^^^")')mix_para*100
!--------------------------------------------------------------------------------------------------------


em_backup=0

!+++++++++++++++++++++SELF CONSISTENCY ITERATION+++++++++++++++++++++++++++++++++++++

temp_E=0.0
open (97,file='../info/rho_error.info')
write(97,'(4X,"ITER NO",10X,"MAX ERROR",10X,"MIN ERROR")')
write(97,'(2X,"==================================================================")')

s_iter=last_iter+1
do iter =s_iter ,max_iter                         ! this is the do loop of iteration

 call cpu_time(start_iter)                   ! Note starting time

 rho_up_old=rho_up		             ! Need to store the old rho s for generation of 	     
 rho_dn_old=rho_dn			     ! new rho s by mixing
 
write(50,*)"STARTING ITERATION  NO :: ",iter
write(*,*)"STARTING ITERATION  NO :: ",iter







!=================CONSTRUCTION OF EFFECTIVE POTENTIAL====================================
!  Vtot_up=Vext+Vh+Vc_up+Vx_up
!  Vtot_dn=Vext+Vh+Vc_dn+Vx_dn
!------------------------------
  call poisson3d(rho_up+rho_dn,vh,E_hart)   	!Genarate hartree POT 
  Vtot_dn=ext_pot+vh 				! Vtot_dn=Vext+Vh 
  Vtot_up=Vtot_dn				! Vtot_up=Vext+Vh
  call add_exchange()  				!Vtot=Vtot+Vx is done in the subroutine
  call add_correlation() 			!Vtot=Vtot+Vc is done in the subroutine
  call write_for_plot (Vtot_up,ndim_z/2,'../output/eff_pot_up.dat')
 call write_for_plot (Vtot_dn,ndim_z/2,'../output/eff_pot_dn.dat')

!========================================================================================






!=============NOW SOLVE SCHRODINGER EQUATION FOR BOTH UP & DN================================
igenE_up=0
igenE_dn=0
Eigen_vec_up=0
Eigen_vec_dn=0
 call solve_3dhamilton(Vtot_up,igenE_up,Eigen_vec_up,1)    
 call solve_3dhamilton(Vtot_dn,igenE_dn,Eigen_vec_dn,2)
!============================================================================================







!===================GENARATE NEW RHO============================================================================================

 call fillup_electrons(occupied_stateup,occupied_statedn)                  ! fillup the electron in the energy levels
 call gen_rho(1,occupied_stateup,Eigen_vec_up,rho_up)			   ! gen rho_up
 call gen_rho(2,occupied_statedn,Eigen_vec_dn,rho_dn)			   ! gen rho_dn
 call check_rho(iter)							   ! check if rho is correctly generated

!===============================================================================================================================






!--------------------------ENERGY CALCULATIONS----------------------------------------------------
 call KE_spin()                                         !calculate Kinetic energy ke_up,ke_dn
 Eext_up=sum(rho_up*ext_pot)				!calculte Eext due to external pot
 Eext_dn=sum(rho_dn*ext_pot)
							!Hartree,exchange & correlation energy E_hart,Ex & Ec are generated before
							!when calculating Vh,Vx & Vc

tot_E=ke_up+ke_dn+Ex+Ec+E_hart+Eext_up+Eext_dn          ! calculate total energy
! call analytic_dft_energy()
!write(50,*)"ENERGY = ",tot_E,"iter =",iter
!---------------------------------------------------------------------------------------------------




!--------A precaution if davidson betrays !!
if((((tot_E-temp_E)/temp_E*100) > max_change_allowed).and. iter>s_iter) then
 write(50,*) " ERROR !!!"
 write(50,*) "Energy in the",iter,"th iteration changed more than the allowed value"
 write(50,*) "It seems that Davidson betrayed us !!!"
 write(50,*) "Doing an energency backup................"
 rho_up=rho_up_old
 rho_dn=rho_dn_old
 call disaster_management(999)
 call charge_den_for_vmd()
 write(50,*) "Program stops here..You can change the davidson's  parameter and continue from this point"
em_backup=em_backup+1
if(em_backup>1)stop
endif



!==========================CONVERGANCE TEST===========================================

! if connvergance achieved then ...............
if (abs(tot_E-temp_E)<=conv_limit) then
  write(50,*)"CONVERGANCE ACHIEVED IN ITERATION NO", iter
  write(50,*) " !!!!!!!!!!!!!   HuRrAy CoNveRgAnCE AchIEvEd  !!!!!!!!!!!!!!!!!!"
  write(50,*)"==================================================================="
  write(6,*) " !!!!!!!!!!!!!   HuRrAy CoNveRgAnCE AchIEvEd  !!!!!!!!!!!!!!!!!!"
  write(6,*)"==================================================================="
  call itime(end_time)
  call cpu_time(sec_end)
 write(60,'(5X,I4,7X,F15.7,5X,F15.7,3X,F15.7,3X,F15.7,5X,F15.7,5X,F10.3)')iter,ke_up+ke_dn,Eext_up+Eext_dn,E_hart,Ex+Ec,tot_E,sec_end-start_iter
  write(60,'(2X,"===============================================================================================================================")')
  write(60,'( "                          !!!!!!!!!!!!!   HuRrAy CoNveRgAnCE AchIEvEd  !!!!!!!!!!!!!!!!!!")')
  write(60,'("                        ========================================================================")')
  write(50,'("PROGRAM ENDS AT ",3I10)')end_time
  write(50,'("TOTALL TIME TAKEN :: ",F20.5," mins")')(sec_end-sec_start)/60.0
  write(60,'("TOTALL TIME TAKEN :: ",F20.5," mins")')(sec_end-sec_start)/60.0
  call write_igenvecup_for_vmd(Eigen_vec_up)          ! these are the subroutines to write the outputs in VMD format
  call write_igenvecdn_for_vmd(Eigen_vec_dn)
  call charge_den_for_vmd()
  call write_for_plot(rho_up,ndim_z/2,'../output/rho_up.dat')
  call write_for_plot(rho_dn,ndim_z/2,'../output/rho_dn.dat')
 stop                          ! It's time to stop
endif




! else we continue ..............

temp_E=tot_E
 call charge_den_mixing(mix_para) 			! Genarate new rho by mixing old & new rho
 call cpu_time(end_iter)				! note end time of iteration
write(60,'(5X,I4,7X,F15.7,5X,F15.7,3X,F15.7,3X,F15.7,5X,F15.7,5X,F10.3)')iter,ke_up+ke_dn,Eext_up+Eext_dn,E_hart,Ex+Ec,tot_E,end_iter-start_iter
 call disaster_management(iter)				! if it's time for backup then do backup else go on ...
  call write_for_plot(rho_up,ndim_z/2,'../output/rho_up.dat')
 call write_for_plot(rho_dn,ndim_z/2,'../output/rho_dn.dat')
 call charge_den_for_vmd()
enddo				! the iteration loop ended atlast !!



!=================================================================================================================================






!if we are here then the convergence is not acheived
  write(50,*)"<<<<<<<<<<<<<<<<<<<<<<<<CONVERGANCE NOT ACHIEVED>>>>>>>>>>>>>>>>>>>>>>>>"
  call write_igenvecup_for_vmd(Eigen_vec_up)
  call write_igenvecdn_for_vmd(Eigen_vec_dn)
  call charge_den_for_vmd()     			!write charge density & spin density in vmd format


end program main


