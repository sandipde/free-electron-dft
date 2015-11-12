subroutine poisson_allocate()
use kind_param
use global_var
use poisson_var
implicit none
!____________ALLOCATIONS_____________________________________________

allocate (d_x(-ndim_x+1:ndim_x))
allocate (d_y(-ndim_y+1:ndim_y))
allocate (d_z(-ndim_z+1:ndim_z))
allocate (double_charge_den(-ndim_x+1:ndim_x,-ndim_y+1:ndim_y,-ndim_z+1:ndim_z))
allocate (gspace_charge_den(ndim_x+1,2*ndim_y,2*ndim_z))
allocate (r_dfg_vhartree  (2*ndim_x, 2*ndim_y, 2*ndim_z))
allocate (c_dfg_vhartree  (ndim_x+1,2*ndim_y,2*ndim_z))
allocate (fft_charge_den(2*ndim_x,2*ndim_y,2*ndim_z))
allocate (tr_x_r2f (-ndim_x+1:ndim_x))
allocate (tr_y_r2f (-ndim_y+1:ndim_y))
allocate (tr_z_r2f (-ndim_z+1:ndim_z))
allocate (tr_x_f2r (2*ndim_x))
allocate (tr_y_f2r (2*ndim_y))
allocate (tr_z_f2r (2*ndim_z))
return
end subroutine poisson_allocate



!_____________________________________________________________________________________________


subroutine poisson_deallocate()
use kind_param
use global_var
use poisson_var
implicit none
!____________DEALLOCATIONS_____________________________________________

deallocate (double_charge_den)
deallocate (gspace_charge_den)
deallocate (r_dfg_vhartree  )
deallocate (c_dfg_vhartree )
deallocate (vhartree)
deallocate (fft_charge_den)
deallocate (tr_x_r2f)
deallocate (tr_y_r2f)
deallocate (tr_z_r2f)
deallocate (tr_x_f2r)
deallocate (tr_y_f2r)
deallocate (tr_z_f2r)
return
end subroutine poisson_deallocate
