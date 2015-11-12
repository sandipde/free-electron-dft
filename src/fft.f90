
!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________
subroutine forward_fft(in,out,dim_x,dim_y,dim_z,delta_r)
use kind_param
!use poisson_var
implicit none
include 'fftw3.f'
integer(I4B)::plan,dim_x,dim_y,dim_z
real(R8B)::in(dim_x,dim_y,dim_z),delta_r
double complex::out(dim_x/2+1,dim_y,dim_z)
out=0.0
!write(*,*)"delta_r=",delta_r
 call dfftw_plan_dft_r2c_3d(plan,dim_x,dim_y,dim_z,in,out,FFTW_ESTIMATE_PATIENT)
 call dfftw_execute(plan)
!write(*,*)plan
 call dfftw_destroy_plan(plan)
out=out*delta_r
return
end subroutine forward_fft

subroutine backward_fft(in,out,dim_x,dim_y,dim_z,cell_vol)
use kind_param
!use poisson_var
implicit none
include 'fftw3.f'
integer(I4B)::invplan,dim_x,dim_y,dim_z
real(R8B)::out(dim_x,dim_y,dim_z),cell_vol,delta_r
double complex::in(dim_x/2+1,dim_y,dim_z)
out=0.0
 call dfftw_plan_dft_c2r_3d(invplan,dim_x,dim_y,dim_z,in,out,FFTW_ESTIMATE)
 call dfftw_execute(invplan)
!write(*,*)invplan
 call dfftw_destroy_plan(invplan)
out=out/cell_vol/8.0
!write(*,*)"cell_vol =",cell_vol
!out=out/(8.0*71*71*71.0)
return
end subroutine backward_fft
