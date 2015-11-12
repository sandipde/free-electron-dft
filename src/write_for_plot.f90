subroutine write_for_plot(func,section,filename)
use kind_param
use global_var
real(R8B)::func(ndim_x,ndim_y,ndim_z)
character(len=50)::filename
integer::i,j,k,section
open(999,file=filename)
do i=1,ndim_x
 do j=1,ndim_y
    write(999,*) xi(i),yi(j),func(i,j,section)
 enddo
 write(999,*)
enddo
close(999)
return
end subroutine write_for_plot

