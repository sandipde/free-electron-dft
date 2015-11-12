subroutine add_myexchange_term()
use kind_param
use global_var
implicit none
real(R8B)::A
real(R8B),allocatable::v_xup(:,:,:),v_xdn(:,:,:),e_x(:,:,:)
allocate (v_xup(ndim_x,ndim_y,ndim_z))
allocate (v_xdn(ndim_x,ndim_y,ndim_z))
allocate (e_x(ndim_x,ndim_y,ndim_z))
A=(3.0/fourpi)**(1.0/3.0)
v_xup=-2.0*A*(rho_up)**(1.0/3.0)
v_xdn=-2.0*A*(rho_dn)**(1.0/3.0)
 call write_for_plot(v_xup,ndim_z/2,"../output/myupexchange_pot.dat")
 call write_for_plot(v_xdn,ndim_z/2,"../output/mydnexchange_pot.dat")
!Vtot_up=Vtot_up-2.0*A*(rho_up)**(1.0/3.0)
!Vtot_dn=Vtot_dn-2.0*A*(rho_dn)**(1.0/3.0)
!allocate (e_x(ndim_x,ndim_y,ndim_z))
e_x=(rho_up**(4.0/3.0)+rho_dn**(4.0/3.0))/(rho_up+rho_dn)
Ex_up=-3.0*A*sum(rho_up*e_x)
Ex_dn=-3.0*A*sum(rho_dn*e_x)
Ex=Ex_up+Ex_dn
deallocate(e_x)
vx=-2.0*A*(rho_up+rho_dn)**(1.0/3.0)
return
end subroutine add_exchange_term