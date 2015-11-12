!_____SANDIP DE_____01.12.08__________________________________________
! This is the matrix vector multiplication routine
! required by dvdson diagonalization subroutine
! This subroutine multiplies the hamilton matrix of original dim N
! (which is stored in compressd row-index storage format )
! with a vector x of dim(N*M1X1) to its right to produce
! C(N*M1,1) vector
!_________________________________________________________________________
 subroutine mat_vec_multi(N,M1,x,C)
use kind_param
use global_var
implicit none
integer(I4B)::M1,N
real(R8B),dimension(N*M1)::x,C
INTEGER(I4B):: i,k,l,index_m1
	C=0.0
do l=1,N*M1
!	index_m1=(l-1)*N
!	do i=1,N  ! n= dim of original hamilton
	!Start with diagonal term.
	i=l-N*((l-1)/(N))
	C(l)=hamilton_elements(i)*x(l)+C(l)
!enddo
!	do l=1,N*M1
!	i=l-N*int(l/(1+N))
	!Loop over oﬀ-diagonal terms.
	index_m1=l-i
	do  k=hamilton_index(i),hamilton_index(i+1)-1
		C(l)=C(l)+hamilton_elements(k)*x(index_m1+hamilton_index(k))
	!	write(*,*)"b(",index(k),")=",b(index(k)),"+,",sa(k)
		C(index_m1+hamilton_index(k))=C(index_m1+hamilton_index(k))+hamilton_elements(k)*x(l)
	enddo 
	enddo

!  do i=1,6
!  write(*,*)C(i)
!  enddo
	return
END subroutine mat_vec_multi
