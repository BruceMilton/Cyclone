      module cyclone_lib
      implicit none
      private
      public :: INTCF,INTCF12,INTCF3,THNRM,PHNRM,allocate_check,AMP,
     +  INTCF3D,INTCFD

! Normalize an angle so that it lies between 0 and 2 PI

      INTERFACE THNRM
         MODULE PROCEDURE THNRM_SINGLE
         MODULE PROCEDURE THNRM_DOUBLE
      END INTERFACE

! Normalize an angle so that it lies between -180 and 180

      INTERFACE PHNRM
         MODULE PROCEDURE PHNRM_SINGLE
         MODULE PROCEDURE PHNRM_DOUBLE
      END INTERFACE


      CONTAINS

!*************************************************************
!
! The subroutine THNRM normalizes an angle in radians so
!  that it lies between 0 and 2 pi
!
      SUBROUTINE THNRM_SINGLE(TH)

      REAL, INTENT(INOUT) :: TH
      REAL, PARAMETER :: TPI=6.28318530717959
      INTEGER I

      do i=1,100
         if (TH >= 0.0) EXIT
         TH = TH +TPI
      end do

      do i=1,100
         if (TH < TPI) EXIT
         TH = TH - TPI
      end do

      return
      END SUBROUTINE THNRM_SINGLE

! ------------------------------------------------------------
      SUBROUTINE THNRM_DOUBLE(TH)

      REAL*8, INTENT(INOUT) :: TH
      REAL*8, PARAMETER :: TPI=6.28318530717959D+0
      INTEGER I

      do i=1,100
         if (TH > 0.0D+0) EXIT
         TH = TH + TPI
      end do

      do i=1,100
         if (TH < TPI) EXIT
         TH = TH - TPI
      end do

      return
      END SUBROUTINE THNRM_DOUBLE
!*************************************************************
!
! PHNRM routine normalizes an angle given in degrees
!  to lie between -180 and +180 degrees
!  An optional second parameter specifies an offset ie.
!   -180 <= (PH -P0) <= 180
!
      SUBROUTINE PHNRM_SINGLE(PH,P0)
      REAL, INTENT(INOUT) :: PH
      REAL, INTENT(IN), OPTIONAL :: P0
      REAL OFFSET

      OFFSET = 0.0
      IF (present(p0)) OFFSET = P0

      DO WHILE ((PH-OFFSET) > 180.0)
         PH = PH - 360.0
      ENDDO

      DO WHILE ((PH-OFFSET) < -180.0)
         PH = PH + 360.0
      ENDDO

      RETURN
      END SUBROUTINE PHNRM_SINGLE
! ------------------------------------------------------------
      SUBROUTINE PHNRM_DOUBLE(PH,P0)
      REAL*8, INTENT(INOUT) :: PH
      REAL*8, INTENT(IN), OPTIONAL :: P0
      REAL*8 OFFSET

      OFFSET = 0.0
      IF (present(p0)) OFFSET = P0

      DO WHILE ((PH-OFFSET) > 180.0)
         PH = PH - 360.0D+0
      ENDDO

      DO WHILE ((PH-OFFSET) < -180.0)
         PH = PH + 360.0D+0
      ENDDO

      RETURN
      END SUBROUTINE PHNRM_DOUBLE
!*************************************************************
!
! Calculate the coefficents for a double 3 point lagrange interpolation
!
      SUBROUTINE INTCF(FR,COF)
      REAL, INTENT(IN) :: FR
      REAL, INTENT(INOUT), DIMENSION(:) :: COF
      REAL :: FR2,FR3

      if (SIZE(COF) < 8) then
         write(*,'('' CALL TO INTCF with vector (COF) size < 8'')')
         stop '*** ILLEGAL CALL TO INTCF ***'
      endif

      FR2=FR*FR
      FR3=FR*FR2
      COF(1) = FR2-.5*(FR3+FR)
      COF(2) = 1.5*FR3-2.5*FR2+1.
      COF(3) = 2.*FR2-1.5*FR3+.5*FR
      COF(4) = .5*(FR3-FR2)
      COF(5) =2.*FR-1.5*FR2-.5
      COF(6) =4.5*FR2-5.*FR
      COF(7) =4.*FR-4.5*FR2+.5
      COF(8) =1.5*FR2-FR
      RETURN
      END SUBROUTINE
!
      SUBROUTINE INTCFD(FR,COF)
      REAL*8, INTENT(IN) :: FR
      REAL*8, INTENT(INOUT), DIMENSION(:) :: COF
      REAL*8 :: FR2,FR3

      if (SIZE(COF) < 8) then
         write(*,'('' CALL TO INTCF with vector (COF) size < 8'')')
         stop '*** ILLEGAL CALL TO INTCF ***'
      endif

      FR2=FR*FR
      FR3=FR*FR2
      COF(1) = FR2-.5*(FR3+FR)
      COF(2) = 1.5*FR3-2.5*FR2+1.
      COF(3) = 2.*FR2-1.5*FR3+.5*FR
      COF(4) = .5*(FR3-FR2)
      COF(5) =2.*FR-1.5*FR2-.5
      COF(6) =4.5*FR2-5.*FR
      COF(7) =4.*FR-4.5*FR2+.5
      COF(8) =1.5*FR2-FR
      RETURN
      END SUBROUTINE INTCFD      
      
      SUBROUTINE INTCF12(FR,COF)
      REAL*8, INTENT(IN) :: FR
      REAL*8, INTENT(INOUT), DIMENSION(:) :: COF
      REAL*8 :: FR2,FR3

      if (SIZE(COF) < 8) then
         write(*,'('' CALL TO INTCF with vector (COF) size < 8'')')
         stop '*** ILLEGAL CALL TO INTCF ***'
      endif

      FR2=FR*FR
      FR3=FR*FR2
      COF(1) = FR2-.5*(FR3+FR)
      COF(2) = 1.5*FR3-2.5*FR2+1.
      COF(3) = 2.*FR2-1.5*FR3+.5*FR
      COF(4) = .5*(FR3-FR2)
      COF(5) =2.*FR-1.5*FR2-.5
      COF(6) =4.5*FR2-5.*FR
      COF(7) =4.*FR-4.5*FR2+.5
      COF(8) =1.5*FR2-FR
      cof(9)=1.0D+0-FR
      cof(10)=-2.0D+0+3.0D+0*FR
      COF(11)=1.0D+0-3.0D+0*FR
      COF(12)=FR
      RETURN
      END SUBROUTINE
!*********************************************************
!
! Calculate the coefficents for a single 3 point lagrange interpolation
!
      SUBROUTINE INTCF3(F,COF)
      real, intent(INOUT) :: F
      real,DIMENSION(:),Intent(INOUT) :: COF

      real :: f2

      if (SIZE(COF) < 7) then
         write(*,'('' CALL TO INTCF3 with vector (COF) size < 7'')')
         stop '*** ILLEGAL CALL TO INTCF3 ***'
      endif

      F2=F*F
      COF(1)=.5*(F2-F)
      COF(2)=1-F2
      COF(3)=.5*(F2+F)
      COF(5)=F-.5
      COF(6)=-2.*F
      COF(7)=F+.5
      RETURN
      END SUBROUTINE INTCF3
!
      SUBROUTINE INTCF3D(F,COF)
      real*8, intent(INOUT) :: F
      real*8,DIMENSION(:),Intent(INOUT) :: COF

      real*8 :: f2

      if (SIZE(COF) < 7) then
         write(*,'('' CALL TO INTCF3 with vector (COF) size < 7'')')
         stop '*** ILLEGAL CALL TO INTCF3 ***'
      endif

      F2=F*F
      COF(1)=.5*(F2-F)
      COF(2)=1-F2
      COF(3)=.5*(F2+F)
      COF(5)=F-.5
      COF(6)=-2.*F
      COF(7)=F+.5
      RETURN
      END SUBROUTINE INTCF3D      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! THERE IS NOW AN INTRINSIC FUNCTION dot_product
!      real*8 function dot(a,b)
!      implicit real*8 (a-h,o-z)
!      real*8 a(3),b(3),dot
!      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
!      return
!      end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc        
! this could be replaced with MATMUL
      subroutine cross(a,b,c)
      implicit real*8 (a-h,o-z)
      real*8 a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine algebra(a,x,b,y,z)
      implicit real*8 (a-h,o-z)
      real*8 x(3),y(3),z(3)
      z(1)=a*x(1)+b*y(1)
      z(2)=a*x(2)+b*y(2)
      z(3)=a*x(3)+b*y(3)
      return
      end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      real*8 function amp(x)
      implicit real*8 (a-h,o-z)
      real*8 x(3)
      amp=sqrt(x(1)**2+x(2)**2+x(3)**2)
      return
      end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smult(scalar,vector_in,vector_out)
!c
!c multiply the input vector by a scalar are return answer
!c  in the output vector
!c
      implicit real*8 (a-h,o-z)
      dimension vector_in(3),vector_out(3)

      vector_out(1)=scalar*vector_in(1)
      vector_out(2)=scalar*vector_in(2)
      vector_out(3)=scalar*vector_in(3)
      return
      end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vnorm(vector)

!c normalize the vector
!c  the function amp is used to find the amplitude of the vector
!c  then each of the components is divided by the amplitude so
!c  that the vector has unit length at return

      implicit real*8 (a-h,o-z)

      dimension vector(3)

      vamp=amp(vector)
      vector(1)=vector(1)/vamp
      vector(2)=vector(2)/vamp
      vector(3)=vector(3)/vamp
      return
      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c  process the allocation stat variable

      subroutine allocate_check(ierr,message)
      integer, intent(IN) :: ierr
      character*(*), intent (IN) ::  message
      
      integer ichar

      if (ierr .ne. 0) then
         write(6,'('' Problem with dynamic allocation '',i5)')ierr
         ichar=len(message)
         write(6,'(a<ichar>)')message
         stop '*** Allocation Error ***'
      endif
      
      return
      end subroutine allocate_check

      end module cyclone_lib
