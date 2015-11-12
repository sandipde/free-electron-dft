subroutine ext_potential()
use kind_param
use global_var
implicit none
integer(I4B)::i,j,k,pot_choice,well_no 
real(R8B)::w(3)       ! wx,wy,wz  for harmonic pot
real(R8B)::xwidth,ywidth,zwidth ! for rectangular well or barrier
real(R8B)::r,rad			! for spherical well
real(R8B)::v0
integer(I4B)::center1,separation			!for double wells
open(2,file='../input/ext_pot.in')
read(2,*)well_no
if (well_no==1)then
read(2,*)pot_choice
eff_vol=cell_vol
!___________Harmonic potential_______________________
if (pot_choice==1)then
	read(2,*) (w(i),i=1,3)
	do i=1,ndim_x
	do j=1,ndim_y
	do k=1,ndim_z
		ext_pot(i,j,k)=0.5*(w(1)*w(1)*xi(i)*xi(i)+w(2)*w(2)*yi(j)*yi(j)+w(3)*w(3)*zi(k)*zi(k))
		enddo
	enddo
	enddo
eff_vol=cell_vol
endif
!__________BOX___________________________________________
if (pot_choice==2) then
	read(2,*) xwidth,ywidth,zwidth
	read(2,*) v0
	ext_pot=v0
	do i=(ndim_x-xwidth)/2,(ndim_x+xwidth)/2
	do j=(ndim_y-ywidth)/2,(ndim_y+ywidth)/2
	do k=(ndim_z-zwidth)/2,(ndim_z+zwidth)/2
		ext_pot(i,j,k)=0.0
	enddo
	enddo
	enddo
eff_vol=(xwidth+1.5)*delta_x*(ywidth+1.5)*delta_y*(zwidth+1.5)*delta_z
write(50,'(3X,"=========================EXTERNAL POTENTIAL============================")')
write(50,'(10X,"TYPE = BOX")')
write(50,'(10X, "DIMENSIONS =",3(2X,F10.4)," Ang")')(xwidth+1.5)*delta_x*bohr2ang,(ywidth+1.5)*delta_y*bohr2ang,(zwidth+1.5)*delta_z*bohr2ang
write(50,'(10X,"VOLUME = ",F10.4,2X,"Ang^3")')eff_vol*(bohr2ang)**3
write(50,'(10X,"POTENTIAL OUTSIDE THE BOX = ",F10.4," Hartree")')v0 
write(50,'(3X,"=======================================================================")')
endif

!______________________SPHERICAL WELL____________________________
if (pot_choice==3) then
	read(2,*)
	read(2,*)v0
	read(2,*)rad
	rad=rad*delta_z
	ext_pot=v0
	do i=1,ndim_x
 		do j=1,ndim_y
			do k=1,ndim_z
			r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
			if (r<=rad)ext_pot(i,j,k)=0.0
			enddo
		enddo
	enddo
eff_vol= 4.0/3.0*pi*rad**3
write(50,'(3X,"=========================EXTERNAL POTENTIAL============================")')
write(50,'(10X,"TYPE = SPHERE")')
write(50,'(10X, "RADIUS =",2X,F10.4," Ang")')rad*bohr2ang
write(50,'(10X,"VOLUME = ",F10.4,2X,"Ang^3")')eff_vol*(bohr2ang)**3
write(50,'(10X,"POTENTIAL OUTSIDE THE SPHERE = ",F10.4," Hartree")')v0 
write(50,'(3X,"=======================================================================")')
endif
!_____________________GAUSSIAN__________________________
if (pot_choice==4)then
	read (2,*)(w(i),i=1,3)
	read(2,*)v0
	do i=1,ndim_x
	do j=1,ndim_y
	do k=1,ndim_z
		ext_pot(i,j,k)=-v0*dexp(-(w(1)*w(1)*xi(i)*xi(i)+w(2)*w(2)*yi(j)*yi(j)+w(3)*w(3)*zi(k)*zi(k)))
		enddo
	enddo
	enddo
endif
!______________________SPHERICAL HARMONIC  WELL____________________________
if (pot_choice==5) then
	read(2,*)(w(i),i=1,3)
	read(2,*)v0
	read(2,*)rad
	rad=rad*delta_z
	ext_pot=0.5*w(1)*w(1)*rad*rad
	do i=1,ndim_x
 		do j=1,ndim_y
			do k=1,ndim_z
			r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
			if (r<=rad)ext_pot(i,j,k)=0.5*w(1)*w(1)*r*r
			enddo
		enddo
	enddo
endif
endif

!===============================DOUBLE WELLS=========================================

if (well_no==2)then
! open(11,file='input/double_pot.in')
read(2,*)pot_choice
!__________BOX___________________________________________
if (pot_choice==2) then
	read(2,*) xwidth,ywidth,zwidth
	read(2,*) v0
!	read(2,*)separation
!	if ((2*xwidth+2*separaton)>=ndim_x)then
 !	   write(*,*)"insufficient points to genarate required potential"
!	   stop
!	endif
	ext_pot=v0
	do i=(ndim_x/2-xwidth)/2,(ndim_x/2+xwidth)/2
	do j=(ndim_y-ywidth)/2,(ndim_y+ywidth)/2
	do k=(ndim_z-zwidth)/2,(ndim_z+zwidth)/2
		ext_pot(i,j,k)=0.0
	enddo
	enddo
	enddo

	do i=(3*ndim_x/2-xwidth)/2,(3*ndim_x/2+xwidth)/2
	do j=(ndim_y-ywidth)/2,(ndim_y+ywidth)/2
	do k=(ndim_z-zwidth)/2,(ndim_z+zwidth)/2
		ext_pot(i,j,k)=0.0
	enddo
	enddo
	enddo
endif

!______________________SPHERICAL WELL____________________________
if (pot_choice==3) then
	read(2,*)
	read(2,*)v0
	read(2,*)rad
	rad=rad*delta_z
!	read(2,*)center1
	read(2,*)separation
	center1=separation/2
	ext_pot=v0
	do i=1,ndim_x
 		do j=1,ndim_y
			do k=1,ndim_z
			r=dsqrt((xi(i))**2+yi(j)*yi(j)+zi(k)*zi(k))
			if (r<=rad) then 
			  ext_pot(i-center1,j,k)=0.0
			  ext_pot(i-center1+separation,j,k)=0.0
			endif
			enddo
		enddo
	enddo
endif

!____________________DAUBLE_GAUSSIAN__________________________
if (pot_choice==4)then
	read (2,*)(w(i),i=1,3)
	read(2,*)v0
	read(2,*)rad	
!	read(2,*)center1
	read(2,*)separation
	center1=separation/2
!	write(*,*)"center = ",center1
	rad=rad*delta_z
	ext_pot=0.0
	do i=1,ndim_x
	do j=1,ndim_y
	do k=1,ndim_z
		if((i-center1)>0)ext_pot(i-center1,j,k)=-v0*dexp(-(w(1)*w(1)*xi(i)*xi(i)+w(2)*w(2)*yi(j)*yi(j)+w(3)*w(3)*zi(k)*zi(k)))+ext_pot(i-center1,j,k)
		if(((i-center1+separation)<=ndim_x).and.((i-center1+separation)>0))then
			ext_pot(i-center1+separation,j,k)=-v0*dexp(-(w(1)*w(1)*xi(i)*xi(i)+w(2)*w(2)*yi(j)*yi(j)+w(3)*w(3)*zi(k)*zi(k)))+ext_pot(i-center1+separation,j,k)
		endif
		enddo
	enddo
	enddo

endif
!______________________SPHERICAL HARMONIC  WELL____________________________
if (pot_choice==5) then
	read(2,*)(w(i),i=1,3)
	read(2,*)v0
	read(2,*)rad
!	read(2,*)center1
	read(2,*)separation
	center1=int(separation)/2
!	center1=int(ndim_x/2-separation/2-rad)
	rad=rad*delta_z
	ext_pot=0.5*w(1)*w(1)*rad*rad
	do i=1,ndim_x
 		do j=1,ndim_y
			do k=1,ndim_z
			r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
			if (r<=rad)then
				ext_pot(i-center1+separation,j,k)=0.0
				ext_pot(i-center1,j,k)=0.0
				ext_pot(i-center1,j,k)=0.5*w(1)*w(1)*r*r+ext_pot(i-center1,j,k)
				ext_pot(i-center1+separation,j,k)=0.5*w(1)*w(1)*r*r+ext_pot(i-center1+separation,j,k)
!			else
!				if ((i-center1)>0)ext_pot(i-center1,j,k)=0.5*w(1)*w(1)*rad*rad
!				if ((i-center1+separation)<=ndim_x)ext_pot(i-center1+separation,j,k)=0.5*w(1)*w(1)*rad*rad
			endif
			enddo
		enddo
	enddo
endif
endif	
1 format (4(F10.3,2X))
open (5,file='../output/pot.dat')
do i=1,ndim_x
	do j=1,ndim_y
		write(5,1) xi(i),yi(j),zi(ndim_z/2),ext_pot(i,j,ndim_z/2)
	enddo
	write(5,*)"   "
enddo
	close(5)
!stop
return
end subroutine ext_potential
