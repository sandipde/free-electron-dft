subroutine calc_rs()
use kind_param
use global_var
implicit none
rs=(3.0/(4.0*pi*n_e/eff_vol))**(1.0/3.0)
write(50,'(15X,"rs = ",F10.4)')rs
return
end subroutine calc_rs
