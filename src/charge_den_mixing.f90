subroutine charge_den_mixing()
use kind_param
use global_var
implicit none
!real(R8B)::mix_para
rho_up=mix_para*rho_up+(1.0-mix_para)*rho_up_old
rho_dn=mix_para*rho_dn+(1.0-mix_para)*rho_dn_old
!write(50,'("^^^^^^^^^^^^MIXING ",F5.2," % OF NEW RHO WITH OLD ONE ^^^^^^^^^^^^^^^^")')mix_para*100
return
end subroutine charge_den_mixing

