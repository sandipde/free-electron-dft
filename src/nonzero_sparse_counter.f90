!________SANDIP DE__01.12.08______________________________________________________
!This subroutine counts the number of nonzero elements of the hamilton 
! the number is required to determine the dimension of the compressed hamilton
! if there are N non zero elements then the hamilton is compressed in 2 one dimensional
! array of length (N+2) 
!___________________________________________________________________________________
subroutine nonzero_sparse_counter()
use kind_param
use global_var
implicit none
integer(i4b)::i,j,k,l,r0,a
allocate (nonzero_counter(ndim_x*ndim_y*ndim_z))
nonzero_counter=1  ! for diagonal elements
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z	
      r0=((i-1)*ndim_y+j-1)*ndim_z+k
      do a=1,max_coeff
	if ((i+a)<=ndim_x) then
		nonzero_counter(r0)=nonzero_counter(r0)+1
        endif
	
        if ((j+a)<=ndim_y)then
		nonzero_counter(r0)=nonzero_counter(r0)+1
	endif
	
        if ((k+a)<=ndim_z)then
		nonzero_counter(r0)=nonzero_counter(r0)+1
	endif
	
      enddo
   enddo
  enddo
enddo
totall_nonzero=sum(nonzero_counter)
deallocate(nonzero_counter)
 compressed_dim=totall_nonzero+2
!write(*,*)"Total Non zero elements:",totall_nonzero
END subroutine nonzero_sparse_counter
