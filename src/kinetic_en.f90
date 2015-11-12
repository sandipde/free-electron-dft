!============================================
! SUBROUTINE TO FIND KE <Psi|del2|Psi>
! SANDIP DE   16.02.09
!=============================================
subroutine kinetic_en(Psi,energy)
use kind_param
use global_var
implicit none
real(r8B)::Psi(Ndim_x*Ndim_y*Ndim_z),energy,intigrand(Ndim_x*Ndim_y*Ndim_z)
 call del2Psi(Ndim_x*Ndim_y*ndim_z,Psi,intigrand)
energy=sum(intigrand*Psi)
return
end subroutine kinetic_en

