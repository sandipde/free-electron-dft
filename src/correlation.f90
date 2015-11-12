subroutine add_correlation()
use kind_param
use global_var
use libxc_funcs_m  
use xc_f90_lib_m
use xc_f90_types_m
implicit none
real(R8B),allocatable::v_cup(:,:,:),v_cdn(:,:,:),e_c(:,:,:)
real(R8B)::temp_rho(2),temp_v(2)
integer(I4B)::i,j,k
TYPE(xc_f90_func_t) :: xc_c_func
TYPE(xc_f90_info_t) :: xc_c_info
allocate (v_cup(ndim_x,ndim_y,ndim_z))
allocate (v_cdn(ndim_x,ndim_y,ndim_z))
allocate (e_c(ndim_x,ndim_y,ndim_z))

 CALL xc_f90_lda_init(xc_c_func, xc_c_info, XC_LDA_C_OB_PW, XC_POLARIZED)
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
     temp_rho(1)=rho_up(i,j,k)
     temp_rho(2)=rho_dn(i,j,k) 
 ! THIS IS SURPRISING THAT ALTHOUGH WE HAVE TO PASS TEMP_RHO(1) ,THE PROGRAM ACTUALLY TAKES TEMP_RHO(2) TOO , IF WE USE SPIN POLARIZED  
      CALL xc_f90_lda_vxc(xc_c_func,temp_rho(1), e_c(i,j,k), temp_v(1))
     v_cup(i,j,k)=temp_v(1)
     v_cdn(i,j,k)=temp_v(2)
    enddo
  enddo
enddo
Vtot_up=Vtot_up+v_cup
Vtot_dn=Vtot_dn+v_cdn
Ec=sum((rho_up+rho_dn)*e_c)     !*delta_x*delta_y*delta_z
!Ec=sum(e_c)
!write (*,*)"Ex=" ,Ex
 call write_for_plot(v_cup,ndim_z/2,"../output/upcorr_pot.dat")
 call write_for_plot(v_cdn,ndim_z/2,"../output/dncorr_pot.dat")
 !call write_for_plot(e_c,ndim_z/2,"../output/e_xc.dat")
deallocate(e_c,v_cup,v_cdn)
return
end subroutine add_correlation
