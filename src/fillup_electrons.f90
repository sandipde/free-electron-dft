subroutine fillup_electrons(occupied_stup,occupied_stdn)
use kind_param
use global_var
implicit none
integer(I4B)::i,j,k,spin,rem_el,occupied_stup,occupied_stdn,degenaracy(nume),ne
real(R8B)::mu,eigval(nume)
fe=0.0
do spin=1,2
  if (spin==1) then
    rem_el=n_eup
    eigval=igenE_up
    ne=n_eup
    degenaracy=0
    call find_degen(eigval,degenaracy)
!    write(*,*)"deg =",degenaracy
    call gen_mu(igenE_up,n_eup,mu)
    
 else
    rem_el=n_edn
    ne=n_edn
    eigval=igenE_dn
    degenaracy=0
    call find_degen(eigval,degenaracy)
    call gen_mu(igenE_dn,n_edn,mu)
    
  endif
 
!write(*,*)"up prob =",fe(1,:)
!write(*,*)"dn prob = ",fe(2,:)
!write(*,*)"degena =",degenaracy
  k=1
 !write(*,*)"rem_el out =",rem_el
 do while (rem_el>0)
    
    if (degenaracy(k) <=rem_el) then
     fe(spin,k)=1.0
!     write (*,*)k,rem_el
     do j=1,degenaracy(k)
      fe(spin,k+j-1)=1.0
!      write (*,*)k+j-1,rem_el
     enddo
!     k=k+1
   else
!     write(*,*)"deg =",degenaracy(k),k
     do j=1,degenaracy(k)
 !    call gen_mu(k,eigval,rem_el,ne,degenaracy,mu)
!     write(*,*)"123"
 !    call fermi_occupancy(eigval(k+j-1),mu,fe(spin,k+j-1))
  !   write(*,*)"fermi =",fe(spin,k+j-1),k
      fe(spin,k+j-1)=rem_el*1.0/degenaracy(k)
   !   write(*,*)"dist =",fe(spin,k+j-1),k
     enddo
!     k=k+degenaracy(k)
   endif
   rem_el=rem_el-degenaracy(k)
    k=k+degenaracy(k)
 enddo
enddo
open (10,file='../info/fermi_occupancy.info')
write(10,'("up prob =",15F15.10)')fe(1,:)
write(10,'("dn prob = ",15F15.10)')fe(2,:)
occupied_stup=0
occupied_stdn=0
do i=1,nume
if (fe(1,i)>=1.0E-5) occupied_stup=occupied_stup+1
if (fe(2,i)>=1.0E-5) occupied_stdn=occupied_stdn+1
enddo
 
write(10,*)"ne=",n_eup,n_edn,"occupied = ",occupied_stateup,occupied_statedn
if (abs(sum(fe(1,:)-n_eup)>0.001)) then
 write(*,*)"fermi occupancy error for up spin"
 stop
endif
if (abs(sum(fe(2,:)-n_edn)>0.001)) then
 write(*,*)"fermi occupancy error for down spin"
 stop
endif
return
end subroutine fillup_electrons

!_______________________________________________________________________________________

subroutine find_degen(energy_spectrum,degenaracy)
use kind_param
use global_var
implicit none
integer(I4B)::i,j,k,degenaracy(NUME),deg
real(R8B)::testE,energy_spectrum(NUME)
degenaracy=0
do i=1,nume
 deg=0
 testE=energy_spectrum(i)
 do j=1,nume
 if (dabs (testE - energy_spectrum(j)) <=deg_limit) deg= deg+1
enddo
degenaracy(i)=deg
enddo
!write(*,*)"degen=",degenaracy
return
end subroutine find_degen

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine fermi_occupancy(E,mu,f)
use kind_param
use global_var
implicit none
real (R8B)::E,f,mu
f=1.0/(1.0+dexp((dabs(E)-mu)/kT))
!write(*,*)"E=",E,"mu=",mu,"F=",f
return
end subroutine fermi_occupancy

!===========================================================================================================================

subroutine gen_mu(eigval,ne,mu)
use kind_param
use global_var
implicit none
real(R8B)::eigval(nume),mu_arr(nume),mu ! the eigen energies
integer(I4B)::state,i,no,n,ne,deg(nume) !n=number of rem electron,ne=tot electron
no=ne
!write(*,*)"deg =",deg
!write(*,*)"no=",no
!state =0
!do state=1,nume
!if (state < n)then
no=state
!if (deg(state)>n)then
! do i=1,deg(state)
!   if ((eigval(state+i)-eigval(state))<=0.001)no=no+1
! enddo
!endif
!endif
!write(*,*)"for mu ne=",no,state
!mu_arr(state)=dabs(eigval(no))
!enddo 
mu=dabs(eigval(ne))
return
end subroutine gen_mu





