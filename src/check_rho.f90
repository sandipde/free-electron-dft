
subroutine check_rho(iter)
use kind_param
use global_var
implicit none
Real(R8B)::intg
integer(I4B)::iter
!****************CHECK IF RHO GENERATION IS CORRECT*********************************
open (1,file='../info/rho_check.info')
 call int_rho(rho_up,intg)
 if (abs(n_eup-intg)<=0.0001)then
    write(1,*)" ALL UP ELECTRONS ARE STILL IN THE SYSTEM....iteration =",iter
 else
   write(1,*) " UP ELECTRONS ARE MISSING/INCREASED !!!"
   write(1,*)"TOTAL UP ELECTRON =",n_eup
   write(1,*)"PRESENT ELECRRON =",intg
   stop
 endif
 call int_rho(rho_dn,intg)
 if (abs(n_edn-intg)<=0.1)then
    write(1,*)" ALL DOWN ELECTRONS ARE STILL IN THE SYSTEM....iteration =",iter
 else
   write(1,*) " DOWN ELECTRONS ARE MISSING/INCREASED !!!"
   write(1,*)"TOTAL DOWN ELECTRON =",n_edn
   write(1,*)"PRESENT ELECRRON =",intg
   write(50,*) "RHO_CHECK ERROR (X)"
   stop
 endif
!************************************************************************************

!--------------------------------------------------------------------------------------
return
end subroutine check_rho

