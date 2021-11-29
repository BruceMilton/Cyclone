! This is a Sample of the file to be used to define a custom ROI routine

Subroutine  Custom_Roi(theta,r,z,tau,nturn,iroi)

use roi
use bfields, only: bcon,bz,dbzdr,dbzdt

implicit none

! don't touch this interface information
! as it is a duplicate of the interface defined in cyclone

real*8,intent(in) :: theta !in radians
real*8,intent(in) :: r
real*8,intent(in) :: z
real*8,intent(in) :: tau ! in radians
integer,intent(in) ::nturn,iroi

!end interface info

! local variables
real :: v
logical :: ifirst=.true.

! General Notes
!  E fields should be in MegaVolts/length unit
!  B fields must be divided by BCON
!  The logical variable INEFLD should be set or cleared for all devices
!   
!  variables available in the ROI module are
!    roi_theta1(nroi) - real*4, entrance angle of roi
!    roi_theta2(nroi) - real*4, exit angle of roi
!    roi_var(nroi) - real*4, user defined value for each roi
!    roi_type(nroi) - integer, type parameter for each Roi, 3=custom
!    nroi - integer, number of defined roi
!    current_roi - integer, this duplictes the value in iroi - don't use
!    inefld - logical, true indicates that the electric field in use
!    er_roi,et_roi,ez_roi - real*8 electric field components

if(ifirst) then
   write(6,'('' ******Warning Default ROI Routine in Use ****'')')
   ifirst=.false.
endif

if (iroi .eq. 1) then ! electrostatic deflector
   V=ROI_VAR(1)/(1E3) !INPUT IN KV/CM, CHANGE TO MV/CM
   ER_ROI=V
   EZ_ROI=0.
   ET_ROI=0.
   INEFLD=.TRUE.
else if(iroi .eq. 2) then !  THIS IS A FIELD FREE REGION
    BZ=0.
    dbzdr=0.
    dbzdt=0.
    INEFLD=.false.
else
   write(6,'('' No Device defined for roi #'',i5)')iroi
END IF

return

end subroutine
