!```````````````````````````````````````````````````````````````````````
! SUBROUTINE FOR CONSTRUCTING & SOLVING 3D HAMILTONIAN
!PROGRAM BY :: SANDIP DE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine solve_3dhamilton(V,eigval,Eigen_vec,sp)
use kind_param
use global_var
implicit none
integer(I4B)::i,j,k,dvdtime1(3),dvdtime2(3),state,sp !sp=1 for up spin ,2 for dn spin
real(R8B)::stime,etime,V(ndim_x,ndim_y,ndim_z),eigval(nume),Eigen_vec(IWRSZ)
external:: mat_vec_multi 
 call cpu_time(stime) 
allocate(hamilton_elements(compressed_dim))

if(sp==0)then
open (11,file='../info/1steigen.info')
write(11,*)"Solving 1st time to construct initial charge density"
else
 open (11,file='../info/eigen.info')
endif
if(sp==1)write(11,*)"SOLVING FOR UP SPINS"
if(sp==2)write(11,*)"SOLVING FOR DOWN SPINS"
 call gen_hamilton(V)
! write(*,*)"333"
! do i=1,10
! 	write(*,*)hamilton_index(i),hamilton_elements(i)
! enddo
 call itime(dvdtime1)
open (2,file='../info/davidson.info')
write(2,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
write(2,*)"DAVIDSON STARTED @ ",dvdtime1,"Hrs"
write(2,*)"Targeting",MBLOCK,"vectors in each iteration................................."
 CALL DVDSON(mat_vec_multi,NMAX,LIM,hamilton_elements(1:hamilton_dim),                   & 
                  ILOW,IHIGH,ISELEC,NIV,MBLOCK,                 & 
                  CRITE,CRITC,CRITR,ORTHO,MAXITER,              & 
                  Eigen_vec,IWRSZ,IWORK,IIWSZ,HIEND,NLOOPS,NMV,IERR)       
 call dvd_err_chk()
 call itime(dvdtime2)
write(2,*)"DAVIDSON COMPLETED @.",dvdtime2,"Hrs"
write(2,*)"Totall time taken::",dvdtime2-dvdtime1
write(2,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
!-----------------------------------------------------------
2000   FORMAT(//8X,'NO',8X,'Eigenvalues',8X,'Eigval Differences',6X,'Residuals',//(8X,I2,F25.15,2F20.10)) 
!WRITE(6,2000) (j,(Eigen_vec(i), i=NUME*hamilton_dim+j,(hamilton_dim+3)*NUME,NUME),j=1,NUME)
!WRITE(3,3000) ((Eigen_vec(i), i=j,j+5), j=1,NUME*hamilton_dim,hamilton_dim)                

write(11,*)"############################################################################"
write(11,1000) nume
WRITE(11,2000) (j,(Eigen_vec(i), i=NUME*hamilton_dim+j,(hamilton_dim+3)*NUME,NUME),j=1,NUME)
WRITE(11,*)"############################################################################"
1000 Format("_________________NO. OF EIGEN VECTORS YOU HAVE :::__",1X,2I)
!do i=1,ndim_x
!	do j=1,ndim_y
!	do k=1,ndim_z
!	write(3,200)(Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+k+state*hamilton_dim),state=0,0)!,Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+ndim_z/2+hamilton_dim), &
!! &	Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+ndim_z/2+2*hamilton_dim),Eigen_vec(((i-1)*ndim_y+j-1)*ndim_z+ndim_z/2+3*hamilton_dim)
!	enddo
!	write(3,*)" "
!	enddo
!	write(3,*) " "
! enddo
 eigval(1:nume) = Eigen_vec(nume*hamilton_dim+1:nume*hamilton_dim + nume)
! call write_igenvec_for_vmd(Eigen_vec)	
 call cpu_time(etime)
write(2,*)"TOTALL TIME TAKEN::",etime-stime ,"sec"
write(2,*)"                  _____________________________________________________________________"
! do i=1,10
! 	write(*,*)Eigen_vec(NUME*original_dim+i)
! enddo 
deallocate(hamilton_elements)
!deallocate(hamilton_index)
return
end subroutine solve_3dhamilton
