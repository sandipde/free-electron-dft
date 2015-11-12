!_____________SANDIP DE_________01.12.08_____________________-
! Subroutine to check the error code produced by dvdson
!______________________________________________________________
       subroutine dvd_err_chk()
       use kind_param
	use global_var
     implicit none
     integer(i4b) :: tag
     
!   if (ierr .eq. 0) then 
!   write(*,*) "Davidson completed normally. This is iteration", NLOOPS
   if (ierr .lt. 0) then  
   write(50,*)" WARNING!....", abs(ierr)," vectors did not converged"
   stop
   elseif (ierr>0)then
   write(50,*)" WARNING!... there is some incosistancy in Davidson. &
		&error code =", ierr
   stop
   endif 


!   ElseIf (INT( MOD(IERR,  2)/1  ) print*,"Oops! N < LIM"
!   ElseIf (INT( MOD(IERR,  4)/2  ) print*,"Oops! LIM < 1"
!   ElseIf (INT( MOD(IERR,  8)/4  ) print*,"Oops! ISELEC(1)<1, and no range specified"
!   ElseIf (INT( MOD(IERR, 16)/8  ) print*,"Oops! IHIGH > N (in range or ISELEC)"
!   ElseIf (INT( MOD(IERR, 32)/16 ) print*,"Oops! IHIGH < ILOW (Invalid range)"
!   ElseIf (INT( MOD(IERR, 64)/32 ) print*,"Oops! NEIG >= LIM (Too many wanted)"
!   ElseIf (INT( MOD(IERR,128)/64 ) print*,"Oops! Probable duplication in ISELEC"
!   ElseIf (INT( MOD(IERR,256)/128) print*,"Oops! NUME >= LIM (max eigen very far)"
!   ElseIf (INT( MOD(IERR,512)/256) print*,"Oops! MBLOCK is out of bounds"
!   ElseIf (INT( MOD(IERR,1024)/512) print*,"Oops! IWRSZ or IIWSZ is not enough"
!   ElseIf (INT( MOD(IERR,2048)/1024)print*,"Oops! Orthogonalization Failed"
!   ElseIf (INT( MOD(IERR,4096)/2048)print*,"Oops! NLOOPS > MAXITER"
!   Else
!     print*, "Program stops here"
!      STOP 
!   endif 


   return 
   end subroutine dvd_err_chk
