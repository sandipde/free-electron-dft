program lxctest

use libxc_funcs_m  
use xc_f90_lib_m
use xc_f90_types_m
!use xc_types
!use libxc

implicit none
real(8) :: rho(2), e_c, v_c1,rho2(2),v_c(2)
TYPE(xc_f90_func_t) :: xc_c_func
TYPE(xc_f90_info_t) :: xc_c_info
rho(1)=1.1
rho(2)=0.5
!v_c=0.0
 CALL xc_f90_lda_init(xc_c_func, xc_c_info,XC_LDA_C_OB_PW, XC_POLARIZED)
 CALL xc_f90_lda_vxc(xc_c_func,rho(1), e_c, v_c(1))
!rho2(1)=1.0
!rho2(2)=0.5
!CALL xc_f90_lda_vxc(xc_c_func,rho2(1), e_c, v_c1)
!CALL xc_f90_lda_end(xc_c_func)
write(*,*)rho
write(*,*)v_c
write(*,*)e_c
!write(*,*)rho2,v_c1
end program lxctest
