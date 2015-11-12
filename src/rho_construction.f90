!================================================================================================================================

subroutine gen_rho(spin,occupied_state,Eigen_vec,rho)
use kind_param
use global_var
implicit none
real(R8B)::rho(ndim_x,ndim_y,ndim_z),psi2(ndim_x,ndim_y,ndim_z)
integer(I4B)::i,j,k,state,occupied_state,spin ! spin=1 for up 2 for dn
real (R8B)::E,Eigen_vec(IWRSZ)
real(R8B)::eigval(nume) ! the eigen energies
integer(I4B)::n !number of electron
rho=0.0
psi2=0.0
!___________________________________________
!rho=f(e)*|Psi|^2
!-------------------------------------------
do state=0,occupied_state-1  ! state=0 is 1st igen vector
do i=1,ndim_x
 do j=1,ndim_y
   do k=1,ndim_z
        !As  Psi(i,j,k)=Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+k+state*hamilton_dim)
	psi2(i,j,k)= Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+k+state*hamilton_dim)*Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+k+state*hamilton_dim) 
!        rho(i,j,k)= ! Psi(0)<>fe(spin,1)
   enddo
  enddo
enddo
!write(*,*)"123"
rho=rho+psi2*fe(spin,state+1)
psi2=0.0
!write(*,*)"345"
enddo
return
end subroutine gen_rho

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine int_rho(rho,intg)
use kind_param
use global_var
implicit none
real(R8B)::rho(ndim_x,ndim_y,ndim_z),intg
integer(I4B)::i,j
intg=0.0
intg=sum(rho)!*(delta_x*delta_y*delta_z)
!write(6,*)"int{rho} = ",intg
return
end subroutine int_rho

subroutine modpsi2(Eigen_vec,val)
use kind_param
use global_var
implicit none
real(R8B)::Psi(ndim_x,ndim_y,ndim_z),val,Eigen_vec(IWRSZ)
integer(I4B)::i,j,k
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
      Psi(i,j,k)=Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+k+1*hamilton_dim)
    enddo
  enddo
enddo
val=sum(Psi*Psi)!*delta_x*delta_y*delta_z
return
end subroutine modpsi2

