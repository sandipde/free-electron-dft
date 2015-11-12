!___SANDIP DE_________01.12.08_____________________
! This is the subroutine to setup the inputs of davidson
! this reads the input from ''dvd.in" file & sets up the 
! parameters  accordingly
! See dvdson subroutine for explanation of this variables
!______________________________________________________
	subroutine dvd_setup()
	use kind_param
	use global_var
      implicit none
      integer :: i     
  open (unit =25,file="../input/dvdson.in")
        Nmax     =  ndim_x*ndim_y*ndim_z
        read(25,*) NUMEMAX
        LIMmax   =  NUMEMAX+35
       IWRSZ    =  2*Nmax*LIMmax + LIMmax*LIMmax +   &     !Dimension of
                   (NUMEMAX+10)*LIMmax + NUMEMAX         ! work array Q
        IIWSZ    =  6*LIMmax + NUMEMAX   
        LIM     = LIMmax   
        read(25,*)  ILOw
        IHIGH   = NUMEMAX                                               
        NUME    =  NUMEMAX                                              
        read(25,*)  niv                                                
        read(25,*)  CRITE
        read(25,*)  CRITC
        read(25,*)  CRITR
        read(25,*)  MBLOCK
        read(25,*)  ORTHO                                         
        MAXITER = MAX( NUME*40, 400 ) 
	allocate (IWORK(IIWSZ))
	allocate (ISELEC(LIM))
      close (unit=25)
!	write(*,*) "YOU HAVE OPTED FOR:",NUME,"EIGEN VALUES & VECTORS"
      return
      end subroutine dvd_setup
