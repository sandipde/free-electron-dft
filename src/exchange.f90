
subroutine add_exchange()
use kind_param
use global_var
use libxc_funcs_m  
use xc_f90_lib_m
use xc_f90_types_m
implicit none
real(R8B),allocatable::v_xup(:,:,:),v_xdn(:,:,:),e_x(:,:,:)
real(R8B)::temp_rho(2),temp_v(2)
integer(I4B)::i,j,k
TYPE(xc_f90_func_t) :: xc_c_func
TYPE(xc_f90_info_t) :: xc_c_info
allocate (v_xup(ndim_x,ndim_y,ndim_z))
allocate (v_xdn(ndim_x,ndim_y,ndim_z))
allocate (e_x(ndim_x,ndim_y,ndim_z))

 CALL xc_f90_lda_init(xc_c_func, xc_c_info, XC_LDA_X, XC_POLARIZED)
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
     temp_rho(1)=rho_up(i,j,k)
     temp_rho(2)=rho_dn(i,j,k) 
 ! THIS IS SURPRISING THAT ALTHOUGH WE HAVE TO PASS TEMP_RHO(1) ,THE PROGRAM ACTUALLY TAKES TEMP_RHO(2) TOO , IF WE USE SPIN POLARIZED  
      CALL xc_f90_lda_vxc(xc_c_func,temp_rho(1), e_x(i,j,k), temp_v(1))
     v_xup(i,j,k)=temp_v(1)
     v_xdn(i,j,k)=temp_v(2)
    enddo
  enddo
enddo
Vtot_up=Vtot_up+v_xup
Vtot_dn=Vtot_dn+v_xdn
Ex=sum((rho_up+rho_dn)*e_x)  !*delta_x*delta_y*delta_z
!Ex=sum(e_x)
!write (*,*)"Ex=" ,Ex
 call write_for_plot(v_xup,ndim_z/2,"../output/upexchange_pot.dat")
 call write_for_plot(v_xdn,ndim_z/2,"../output/dnexchange_pot.dat")
! call write_for_plot(e_x,ndim_z/2,"../output/e_xc.dat")
deallocate(e_x,v_xup,v_xdn)
return
end subroutine add_exchange
