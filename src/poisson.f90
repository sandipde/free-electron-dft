!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________

subroutine poisson3d(charge_den,hartree_pot,hartree_energy)
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k,ix,iy,iz,i1,choice,no
real(r8B)::exact_energy,rmax,r,max_diff,min_diff,E
real(r8B),allocatable::exact_pot(:,:,:),pot_diff(:)
real(R8B)::charge_den(ndim_x,ndim_y,ndim_z),hartree_pot(ndim_x,ndim_y,ndim_z)
real(R8B)::hartree_energy
!allocate(charge_den(ndim_x,ndim_y,ndim_z))
!allocate (hartree_pot(ndim_x,ndim_y,ndim_z))
! call read_input()


charge_den=charge_den/(delta_x*delta_y*delta_z)         ! The reason for this
!                               is that this poisson solver was written for
!                               unnormalized rho but in our program we have rho
!                               normalized....So without modifying this program
!                               I simply add this line


open (3,file='../info/poisson.info')
 write(3,*)"MIN RAM REQUIRED > ",(3.0*(ndim_x+ndim_y+ndim_z)+35.0*ndim_x*ndim_y*ndim_z+16.0*ndim_y*ndim_z)*8.0/(1024*1024),"MB"
! include 'allocate'
 call poisson_allocate()
open (55,file="../output/charge_den.dat")
do i=1,ndim_x
 do j=1,ndim_y
   write(55,*)xi(i),yi(j),charge_den(i,j,ndim_z/2)
 enddo
   write(55,*)""
enddo
!write(*,*)"alll"
alpha_conv_L=5.0d0
! call gen_mesh()
!write(6,*)"plz enter your choice:: 1= gaussian  2=charged sphere"
! read(*,*) choice
!  if (choice==1)then
 ! call charge_density()
!  exact_energy=dsqrt(2.0d0/pi)*kappa
!  else
!    rmax=min(Lx,Ly,Lz)/2.0
!    write(*,*)"radius=" ,rmax
!   call gen_charged_sphere(rmax)
!   exact_energy=1.0/(pi)/(4.0/3.0*pi*rmax**3)
!  endif
 call double_grid_trans(charge_den)
  
 call fft_grid_trans()
!write(*,*)"123"
 call gen_hartree_pot(hartree_pot,hartree_energy)
 call poisson_deallocate()
write(3,*)"___________________________________OUTPUT_____________________________________"
write(3,'(a35,2x,f15.7)')"Integral rho(r) dr  = ", sum(charge_den)*delta_r
!write(6,'(a35,2x,f15.13)')"Analytical energy ~", exact_energy
write(3,'(a35,2x,f15.10)')"Numerical G-space energy =",hartree_energy !/delta_r
!hartree_energy =sum(charge_den*hartree_pot) * delta_r
!hartree_energy=sum(charge_den*hartree_pot)
write(3,'(a35,2x,f15.10)')"Numerical R-space energy  = ",hartree_energy
!write(3,'(a35,2x,f15.7)')"Error(%) =",(exact_energy - hartree_energy)/exact_energy*100.0d0
write(3,*)"____________________________________END_________________________________________"
! allocate(exact_pot(ndim_x,ndim_y,ndim_z))
!   call analytic_pot(rmax,exact_pot)
!open (5,file="hartree_pot.dat")
! if (choice==1)then
! no=0
! kappa  =  a0/2.0
! allocate (pot_diff(0:ndim_x*ndim_y*ndim_z))
! pot_diff=0.0
! min_diff=0
! max_diff=0
!   do i=1,ndim_x
!     do j=1,ndim_y
!       do k=1,ndim_z
!         no=no+1
!         r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
!         if (r==0.0) then
!           exact_pot(i,j,k)= 2.0d0/dsqrt(pi) * kappa
!         else
!         exact_pot(i,j,k)=derf(kappa*r)/r
!         endif
!        
!        pot_diff(no)=exact_pot(i,j,k)-hartree_pot(i,j,k)
!       enddo
!     enddo
!  enddo
! write(*,*) pot_diff(maxloc(pot_diff)),pot_diff(minloc(pot_diff))
open(500,file='../output/hartree_pot.dat')
100 Format (3F20.15,2x)
  do i=1,ndim_x
  do j=1,ndim_y
       write(500,100) xi(i),yi(j),hartree_pot(i,j,ndim_z/2)!,exact_pot(i,j,ndim_z/2)
   enddo
   write(500,*)" "
  enddo
 close(500)
! else
!   do i=1,ndim_x
!   do j=1,ndim_y
!    do k=1,ndim_z
!        write(5,*) dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)),hartree_pot(i,j,k),hartree_pot(i,j,k)/exact_pot(i,j,k)
!    enddo
!   enddo
!   enddo
! endif
  
  

 !call write_in_vmd_format(hartree_pot,'hartree_vmd.dat')
end subroutine poisson3d
