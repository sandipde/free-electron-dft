!___________________________________________________________________________________________________________________________
! SUBROUTINE TO GENARATE HARTREE POTENTIAL & HARTREE ENERGY IN FOURIER SPACE
!                  PROGRAM BY:: SANDIP DE 26.01.09
!____________________________________________________________________________________________________________________________
subroutine gen_hartree_pot(hartree_pot,hartree_energy)
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k,ix,iy,iz,i1
real(R8B)::z2,yz2,r_sqr,rvalue,error_fn,hartree_pot(ndim_x,ndim_y,ndim_z),hartree_energy
!real(R8B)::hartree_pot(ndim_x,ndim_y,ndim_z),hartree_energy
delta_r=delta_x*delta_y*delta_z
alpha_x = alpha_conv_L/Lx
alpha_y =alpha_conv_L/Ly
alpha_z = alpha_conv_L/Lz
alpha=dsqrt((alpha_x*alpha_x+alpha_y*alpha_y+alpha_z*alpha_z)/3.0d0)
alpha_sqr=alpha*alpha
!write(*,*)"alpha_x=",alpha_x,"alpha=",alpha
!write(*,*)"cell vol=",cell_vol

!_________________________________________________________________________
! first take the fourier transform of the density which is
! stored in "fft_charge_den" in fft mesh...
!__________________________________________________________________________
 call forward_fft(fft_charge_den,gspace_charge_den,2*ndim_x,2*ndim_y,2*ndim_z,delta_r)
 !call backward_fft(gspace_charge_den,fft_charge_den,2*ndim_x,2*ndim_y,2*ndim_z,cell_vol)

!____________________________________________________________________________________________________
! calculate long range part of the 1/r term i.e. Vc_long,
! given as
!
!               erf(\alpha r)
!  Vc_long(r) = -------------
!                    r
! this is done directly on the fft mesh...
! the result is then stored in the array "vhartree"...

!____________________________________________________________________________________________________

do ix =-ndim_x+1,ndim_x
  i      =  tr_x_r2f(ix) 
  do iy= -ndim_y+1,ndim_y
   j   =  tr_y_r2f(iy)
    do iz=-ndim_z+1,ndim_z
       k  =  tr_z_r2f(iz)
       r_sqr=d_x(ix)*d_x(ix)+d_y(iy)*d_y(iy)+d_z(iz)*d_z(iz)
       rvalue = dsqrt(r_sqr)
       error_fn=derf(alpha*rvalue)
       r_dfg_vhartree (i, j, k) = error_fn/rvalue
     enddo
   enddo
enddo
!write(*,*)"qww"
deallocate (d_x)
deallocate (d_y)
deallocate (d_z)
!write(*,*)"rst"
! the r --> zero limit of erf(\alpha r)/r = 2/root(pi) * \alpha
! set "r_dfg_vhartree" at r=0 to this value
!
      i0    =  tr_x_r2f(0)
      j0    =  tr_y_r2f(0)
      k0    =  tr_z_r2f(0)
!write(*,*)"i0,j0,k0=",i0,j0,k0
      r_dfg_vhartree (i0,j0,k0) = 2.0d0/dsqrt(pi) * alpha

! now take the fourier transform of the long range term of 1/r,
! i.e. "Vc_long" computed above...
!
 call forward_fft(r_dfg_vhartree,c_dfg_vhartree,2*ndim_x,2*ndim_y,2*ndim_z,delta_r)


! compute the short range part of the 1/r
! which is analytically available in G-space...
!                       4 \pi               -g^2
! i.e. Vc_short (g)  =  ----- ( 1 - exp (----------)  )
!                        g^2              4\alpha^2
!
! further, since we already have the fourier transformed
! density i.e. \rho(g), we then directly multipy it
! to the sum of the short and long range of 1/r,
! i.e. (Vc_short+Vc_long) which is stored in "vhartree".
! this give the final "vhartree" in G-space and fft mesh.
 call gen_gsqr()
 four_alpha_sqr=4.0d0*alpha*alpha
 vsumg  = 0.0d0
 vtermg = 0.0d0

      do k = 1, 2*ndim_z
         do j = 1, 2*ndim_y
            do i = 1, ndim_x+1
               g_sqr  = gsqr(i,j,k)
               if ( g_sqr .ne. 0.0d0) then
                  g_exp_term = dexp(-g_sqr / four_alpha_sqr)
                  fourpi_by_gsqr = fourpi/g_sqr

! add the short range part of the potential to vhartree(g) = erf/r
                  c_dfg_vhartree(i,j,k) = c_dfg_vhartree(i,j,k)+ fourpi_by_gsqr*(1.0d0-g_exp_term)
! multiply "gspace_charge_den(g)" to vhartree(g) to get the hartree potential
                  c_dfg_vhartree(i,j,k) = c_dfg_vhartree(i,j,k)*(gspace_charge_den(i,j,k))

! multiply conjg of "gspace_charge_den(g)" to the above hartree potential
! vhartree(g) to calculate the hartree energy...
                  vtermg =dconjg(gspace_charge_den(i,j,k))* c_dfg_vhartree(i,j,k)
                  vsumg  = vsumg + vtermg

               endif
            enddo
         enddo
      enddo
!     write(*,*)"vsum=",vsumg
      vtermg = 0.0d0
      do k   = 1, 2*ndim_z
         do j  =  1, 2*ndim_y
            do i  =  ndim_x+2,2*ndim_x
               ix =  tr_x_f2r(i)
               i1 =  tr_x_r2f(-ix)
		!write(*,*)"i=",i,"ix=",ix,"i1=",i1
               vtermg =(gspace_charge_den(i1,j,k))*dconjg (c_dfg_vhartree(i1,j,k))
               vsumg  = vsumg + vtermg
            enddo
         enddo 
      enddo
!  write(*,*)"vsum=",vsumg
!64
! the g --> zero limit of long range part of "1/r" potential
! is -pi/\alpha^2,

      deallocate(gsqr)
      vtermg = 0.0d0
      c_dfg_vhartree (i0,j0,k0) = c_dfg_vhartree (i0,j0,k0)-(-pi/alpha_sqr)
     c_dfg_vhartree (i0,j0,k0) = c_dfg_vhartree (i0,j0,k0)* gspace_charge_den(i0,j0,k0)
      vtermg =dconjg( gspace_charge_den(i0,j0,k0)) * c_dfg_vhartree(i0,j0,k0)
!	write(*,*)"vtermg=",vtermg
      hartree_energy = dreal(vsumg+vtermg) / (8.0d0*cell_vol)
!__________________________________________________________________________________________
! finally, take the inverse fourier transform of the hartree
! potential to get the potential in Real-space_FFT-mesh...
!

      r_dfg_vhartree = 0.0d0
      call backward_fft(c_dfg_vhartree, r_dfg_vhartree,2*ndim_x, 2*ndim_y,2*ndim_z,cell_vol)
allocate (vhartree (-ndim_x/2+1:ndim_x/2, -ndim_y/2+1:ndim_y/2, -ndim_z/2+1:ndim_z/2))
      call double_fg_2_single()
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
      hartree_pot(i,j,k)=vhartree(i-ndim_x/2,j-ndim_y/2,k-ndim_z/2)
    enddo
  enddo
enddo
end subroutine gen_hartree_pot
