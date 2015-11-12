!_____SANDIP DE_____01.12.08__________________________________________
! This is the matrix vector multiplication routine
! This subroutine multiplies the Ke_hamilton matrix of original dim N
! (which is stored in compressd row-index storage format )
! with a vector x of dim(N) to its right to produce
! C(N,1) vector
!_________________________________________________________________________
 subroutine del2Psi(N,x,C)
use kind_param
use global_var
implicit none
integer(I4B)::M1,N
real(R8B),dimension(N)::x,C
INTEGER(I4B):: i,k,l,index_m1
M1=1
do l=1,M1
	index_m1=(l-1)*N
	do i=1,N  ! n= dim of original hamilton
	!Start with diagonal term.
	C(index_m1+i)=ke_hamilton_elements(i)*x(index_m1+i)
	enddo
	do i=1,N
	!Loop over oﬀ-diagonal terms.
	do  k=hamilton_index(i),hamilton_index(i+1)-1
		C(index_m1+i)=C(index_m1+i)+ke_hamilton_elements(k)*x(index_m1+hamilton_index(k))
	!	write(*,*)"b(",index(k),")=",b(index(k)),"+,",sa(k)
		C(index_m1+hamilton_index(k))=C(index_m1+hamilton_index(k))+ke_hamilton_elements(k)*x(index_m1+i)
	enddo 
	enddo
enddo
!  do i=1,6
!  write(*,*)C(i)
!  enddo
	return
END subroutine del2Psi
