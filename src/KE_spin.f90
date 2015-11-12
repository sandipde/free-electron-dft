subroutine KE_spin()
use kind_param
use global_var
implicit none
integer(I4B)::i
real(R8B)::ke
ke_up=0.0
ke_dn=0.0
ke=0.0
do i=0,occupied_stateup-1
 call kinetic_en(Eigen_vec_up(i*hamilton_dim+1:(i+1)*hamilton_dim),ke)
ke_up=ke_up+fe(1,i+1)*ke
enddo
ke=0.0
do i=0,occupied_statedn-1
 call kinetic_en(Eigen_vec_dn(i*hamilton_dim+1:(i+1)*hamilton_dim),ke)
ke_dn=ke_dn+fe(2,i+1)*ke
enddo
write(50,*)"kinetic energy calculated"
write(50,'("KE_up = ",F20.10," KE_dn = ",F20.10)')ke_up,ke_dn
return
end subroutine KE_spin
