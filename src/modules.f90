module kind_param

! * * * * * * * * SYMBOLIC NAMES OF VARIABLE TYPES * * * * * * * *

! 4-, 8-, 2-, and 1-byte integers
  INTEGER, PARAMETER :: I4B=SELECTED_INT_KIND(8)
  INTEGER, PARAMETER :: I8B=SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: I2B=SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B=SELECTED_INT_KIND(2)
  
! 8- and 4-byte reals
  INTEGER, PARAMETER :: R8B=SELECTED_REAL_KIND(10)
  INTEGER, PARAMETER :: R4B=SELECTED_REAL_KIND(4)
  
end module kind_param

module global_var
use kind_param
 character(len=3)::ans             !if ans=y create new density else read density
 character(len=3)::fstop,cont,start_from_backup          !if fstop=y stop after creating the
                                         !density
!__________________SELF CONSIST_______________________________

integer(I4B)::n_e,n_eup,n_edn,backup_iter,last_iter,s_iter
real(R8B),allocatable,dimension(:,:,:)::rho_up_old,rho_dn_old,exact_rho
integer(I4B)::occupied_stateup,occupied_statedn
real(R8B)::ke_up,ke_dn
real(R8B)::mix_para  ! mixing parameter for charge density
real(R8B)::conv_limit,max_change_allowed !convergance energy limit
integer::max_iter
!************MESH*********************************************************

integer(I8B)::ndim_x,ndim_y,ndim_z		                ! grid dimension
real(R8B),DIMENSION(:),allocatable::xi,yi,zi               ! x,y,z points
real(R8B)::delta_x,delta_y,delta_z                                ! Mesh widths
real(R8B)::x_bound,y_bound,z_bound                          ! Bounderies of the Mesh
real(R8B)::Lx,Ly,Lz,cell_vol,rs,eff_vol

!****************************************************************************

!*************Potential*******************************************************
real(R8B),DIMENSION(:,:,:),allocatable::ext_pot,Vtot_up,Vtot_dn,Vtot
!****************************************************************************

!***********Hamiltoninan****************************************************
integer(I4B),allocatable::nonzero_counter(:),hamilton_index(:)              !hamilton_index stores the indexes
integer(I4B)::compressed_dim,totall_nonzero,hamilton_dim,max_coeff
real(R8B),allocatable,dimension(:)::ke_hamilton_elements,hamilton_elements,igenE_up,igenE_dn                                           !stores the elements 

!_______________________________________________________________________________

!------------------------------------------
! DVDSON
   INTEGER(I4B) :: NMAX, NUMEMAX,LIMMAX, IWRSZ,NUME,   &   
                   IIWSZ ,LIM,ILOW,IHIGH, NIV,         &
                   MBLOCK,MAXITER,NLOOPS,NMV,IERR
   REAL(R8B)    :: CRITE , CRITC, CRITR, ORTHO  

   INTEGER(I4B),ALLOCATABLE :: ISELEC(:), IWORK(:)
REAL(R8B), ALLOCATABLE         :: Eigen_vec_up(:),Eigen_vec_dn(:)
   logical ::HIEND

!_____________________________________________________________________
! POISSON
!___________________________________________________________________

real(R8B),allocatable::rho_up(:,:,:),vh_up(:,:,:),rho_dn(:,:,:),vh_dn(:,:,:),vh(:,:,:)
real(R8B)::hart2ev, ev2hart, bohr2ang, ang2bohr, a0,E_hart,E_hart_up,E_hart_dn
parameter  (hart2ev=27.2116d0,   ev2hart=1.0d0/hart2ev)
parameter  (bohr2ang=0.529177d0, ang2bohr= 1.0d0/bohr2ang)
parameter  (a0=0.529177d0)

!========================================================================
!EXCHANGE-CORRELATION POTENTIAL
!-------------------------------------------------------------------------
real(R8B)::Ex_up,Ex_dn,Eext_up,Eext_dn,Ex,Ec
real(R8B),allocatable::vx(:,:,:)
!---------------------------------------------------------------------------
real(R8B)::analytic_E

!______________________________________________________________________________
!   CHARGE DENSITY CONSTRUCTION
!_____________________________________________________________________________
real(R8B),allocatable::fe(:,:) !fermi occupancy fe(1/2,:) = accupancy for spin up/dn
real(R8B)::kT,pi,twopi,fourpi,deg_limit
parameter (pi=3.141592653589793238462643383d0)
parameter(twopi=2.0d0*pi)
parameter(fourpi=4.0d0*pi)	
end module global_var

!################################################################################################################################

module poisson_var
use kind_param
!__________double grid__________________
real(R8B),DIMENSION(:),allocatable::d_x,d_y,d_z
real(R8B)::dx_bound,dy_bound,dz_bound
!_____________Charge densisity_______________________________

real(R8B)::gauss_const,kappa
real(R8B),allocatable::fft_charge_den(:,:,:)
real(R8B),allocatable::double_charge_den(:,:,:)      ! charge density in the double grid
integer(I4B)::original_x_start,original_y_start,original_z_start
!_______fft__________________________________
real(R8B),allocatable,dimension(:)::fft_gx,fft_gy,fft_gz
integer(i4B),allocatable,dimension(:)::tr_x_r2f,tr_y_r2f,tr_z_r2f,tr_x_f2r,tr_y_f2r,tr_z_f2r     !fourier<>real mesh maps
double complex,allocatable::gspace_charge_den(:,:,:),c_dfg_vhartree(:,:,:)
real(R8B),allocatable::r_dfg_vhartree(:,:,:),vhartree (:,:,:)
integer(I4B)::plan,invplan,i0,j0,k0
real(R8B)::delta_r
real(R8B),allocatable:: gsqr (:,:,:)
real(R8B)::tpibyLx, tpibyLy, tpibyLz
real(R8B)::gz2, gyz2, g_sqr, g_sqr_max
real(R8B)::four_alpha_sqr, g_exp_term ,fourpi_by_gsqr,alpha_sqr
!real(R8B)::hartree_energy
double complex:: vtermg, vsumg
real(R8B)alpha_x,alpha_y,alpha_z,alpha,alpha_conv_L
end module poisson_var

!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________
