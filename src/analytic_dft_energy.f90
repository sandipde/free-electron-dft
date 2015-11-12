subroutine analytic_dft_energy()
use kind_param
use global_var
implicit none
real(R8B)::eigen_sumup,eigen_sumdn
integer(I4B)::i,j,k,spin
eigen_sumup=0.0
do i=1,occupied_stateup
 eigen_sumup=eigen_sumup+fe(1,i)*igenE_up(i)
enddo
eigen_sumdn=0.0
do i=1,occupied_statedn
 eigen_sumdn=eigen_sumdn+fe(2,i)*igenE_dn(i)
enddo
analytic_E=eigen_sumup+eigen_sumdn-E_hart+(Ex_up+Ex_dn)-sum(vx*(rho_up+rho_dn))
!write(*,*)eigen_sumup,eigen_sumdn
return
end subroutine analytic_dft_energy
