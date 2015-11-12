!______SANDIP DE_____16.02.09_______________________________________________________
! Subroutine to generate the KE _hamilton in compressed storage format
! The hamilton is symmetric sparse
! So I stored only the upper triangle in row-index storage format
! In this format we only stores the nonzero elements.If N is the totall non-zero
! elements than we store it in two 1D array each of length N+1
! As here I'm storing only the upper triangle I need  two 1D array each of length U+2
! U where is the totall nonzero elements in upper triangle .
!______________________________________________________________________________________
subroutine gen_kinetic_hamilton()
use kind_param
use global_var
implicit none
real(r8b)::inv_delta_x2,inv_delta_y2,inv_delta_z2
integer(i4b)::i,j,k,l,r0,r,a,coeff,stime1(3),etime1(3)
integer(i4b)::counter,index_count,index_off_first
!____cx(n)=coeff of Si((i+n),j,k) in the expansion of second derivative of Si(i,j,k)___
real(r8b)::cx(6),cy(6),cz(6) !expansion coeff for x,y,z part of wave functon
real(r8b)::c0                        ! coeff for diagonal elements
 call itime(stime1)
include "hamilton_coeff"  ! coefficients are evaluated in this file
hamilton_dim=ndim_x*ndim_y*ndim_z
 counter=0
!_______Storing diagonals____________________________________

do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z	
      r0=((i-1)*ndim_y+j-1)*ndim_z+k
	counter=counter+1
 !       write(*,*)"111",i,j,k
	ke_hamilton_elements(counter)=c0     
     enddo
  enddo
enddo
!write(*,*)"qrrrrr"
!____________________________________________________________

index_count=0
index_off_first=0         ! this is index for first N places of hamilton_index
 counter=hamilton_dim+1
ke_hamilton_elements(counter)=0  ! ARBITRARY(according to this storage rule)
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z	
      r0=((i-1)*ndim_y+j-1)*ndim_z+k
	index_off_first=index_off_first+1

	hamilton_index(index_off_first)=counter+1  ! first N places contains the location of first off diagonal
										! elements in hamilton_elements.Here N=hamilton_dim
	 do a=1,max_coeff
	 if ((k+a)<=ndim_z)then
		r=r0+a
		counter=counter+1
		index_count=index_count+1
		ke_hamilton_elements(counter)=cz(a)
		hamilton_index(hamilton_dim+1+index_count)=r
	endif
	enddo
	 do a=1,max_coeff
	 if ((j+a)<=ndim_y)then
		r=r0+ndim_z*a
		counter=counter+1
		index_count=index_count+1
		ke_hamilton_elements(counter)=cy(a)
		hamilton_index(hamilton_dim+1+index_count)=r
	endif
	enddo
	 do a=1,max_coeff
	if ((i+a)<=ndim_x) then
		r=r0+ndim_y*ndim_z*a
		counter=counter+1
		index_count=index_count+1
		ke_hamilton_elements(counter)=cx(a)
		hamilton_index(hamilton_dim+1+index_count)=r
        endif
      enddo
ke_hamilton_elements(compressed_dim)=0					! Taking
hamilton_index(hamilton_dim+1)=compressed_dim+1		! Care of the Rules
hamilton_index(compressed_dim)=hamilton_dim-1			! as I stored only half matrix
   enddo
  enddo
enddo
call itime(etime1)
!write(*,*)"##################################################"
!write(*,*)"                          sparse hamilton  compressed succesfully !!!!!!!"
!write(*,*)"                     original dimension::",hamilton_dim,"X",hamilton_dim
!write(*,*)"               compressed in two 1D arrays of dimesion::",compressed_dim
!write(*,*)"                       compression quality(%)::",float(hamilton_dim*hamilton_dim-compressed_dim*2)/(hamilton_dim*hamilton_dim)*100.0
!write(*,*)"time taken::",etime1-stime1
!write(*,*)"__________________________________________________________________________________"

! open (10,file="hamilton.dat")
! do i=1,compressed_dim
! 	write(10,*) i,hamilton_index(i),hamilton_elements(i)
! enddo

return
END subroutine gen_kinetic_hamilton
