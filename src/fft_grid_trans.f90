!____________________________________________________________________________________
!	PROGRAM BY SANDIP DE   25.01.09
!____________________________________________________________________________________
subroutine fft_grid_trans()
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k

!______________________________________________________________________________________________________________________

!allocate (fft_gx(2*ndim_x))
!allocate (fft_gy(2*ndim_y))
!allocate (fft_gz(2*ndim_z))

!______________________________________________________________________________________________________________________

 call gen_fr_maps(tr_x_r2f,tr_x_f2r,2*ndim_x)
 call gen_fr_maps(tr_y_r2f,tr_y_f2r,2*ndim_y)
 call gen_fr_maps(tr_z_r2f,tr_z_f2r,2*ndim_z)
!write(*,*)"123"
fft_charge_den=0.0
do i=1,2*ndim_x
  do j=1,2*ndim_y
    do k=1,2*ndim_z
	!write(*,*)i,j,k,tr_x_f2r(i),tr_y_f2r(j),tr_z_f2r(k)
      fft_charge_den(i,j,k)=double_charge_den(tr_x_f2r(i),tr_y_f2r(j),tr_z_f2r(k))
    enddo
  enddo
enddo
!write(*,*)"345"
!open (3,file='fft_charge_den.dat')
!k=ndim_z
!!do i=1,2*ndim_x
!  do j=1,2*ndim_y
	!write(*,*)i,j,k,tr_x_r2f(i),tr_y_r2f(j),tr_z_r2f(k)
!    write(3,200)d_x(tr_x_f2r(i)),d_y(tr_y_f2r(j)),d_z(ndim_z),fft_charge_den(i,j,1)
 ! enddo 
  !  write(3,*)"  "
!enddo
!200 Format (4(F15.7,2X))
return
end subroutine fft_grid_trans

!__________________________________________________________________________________________________________

! SUBROUTINE TO GEN MAPS BETWEEN REAL AND FOURIER MESH

!__________________________________________________________________________________________________________
 subroutine gen_fr_maps (tr_r2f, tr_f2r, tr_dim)
 use kind_param
 implicit none
      integer (I4B)::  tr_dim
      integer(I4B):: tr_r2f(-tr_dim/2+1:tr_dim/2),tr_f2r(tr_dim)
      integer (I4B)::  i, j, tr_dim_by2m1, tr_dim_by2

!--------------------------
  
      tr_dim_by2   = tr_dim/2
      tr_dim_by2m1 = tr_dim_by2-1

      tr_r2f  =  0
      tr_f2r  =  0
	j=0
      do i = -tr_dim_by2m1, tr_dim_by2
         if(i.ge.0) then
            j = i + 1
         else
            j = i + tr_dim + 1
         endif
         tr_r2f(i)  =  j
         tr_f2r(j)  =  i
      enddo
!write(*,*)"in",tr_r2f
    return
      end subroutine gen_fr_maps

!______________________________________________________________________________________________________________
!  SUBROUTINE TO GENERATE gsqr
!_______________________________________________________________________________________________________________

subroutine gen_gsqr()
use kind_param
use global_var
use poisson_var
implicit none
integer(i4B)::i,j,k,ix,iy,iz
!--------------------------
allocate (gsqr (2*ndim_x, 2*ndim_y,2*ndim_z))
      gsqr    = 0.0d0
      tpibyLx = twopi/(2*Lx)
      tpibyLy = twopi/(2*Ly)
      tpibyLz = twopi/(2*Lz)

!      g_sqr_max = 0.0d0
      do k   = 1, 2*ndim_z
         iz  = tr_z_f2r(k)
         gz2 = (iz*tpibyLz)*(iz*tpibyLz)
         do j    = 1, 2*ndim_y
            iy   = tr_y_f2r(j)
            gyz2 = (iy*tpibyLy)*(iy*tpibyLy) + gz2
            do i     = 1, 2*ndim_x
               ix    = tr_x_f2r(i)
               gsqr(i,j,k) = (ix*tpibyLx)*(ix*tpibyLx) + gyz2
            enddo
         enddo
      enddo
return
end subroutine gen_gsqr

subroutine double_fg_2_single()
use kind_param
use poisson_var
use global_var
implicit none
integer(I4B)::i,j,k,ix,iy,iz
vhartree=0.0
do ix=-ndim_x/2+1,ndim_x/2
  i  = tr_x_r2f(ix)
  do iy=-ndim_y/2+1,ndim_y/2
  j = tr_y_r2f(iy)
    do iz=-ndim_z/2+1,ndim_z/2
    k=tr_z_r2f(iz)
    vhartree(ix,iy,iz)=r_dfg_vhartree(i,j,k)
    enddo
  enddo
enddo
return
end subroutine double_fg_2_single


















