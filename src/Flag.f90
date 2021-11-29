module Flag
      use cyclone_lib
      use cyclone_data
      use io
      implicit none      
      private

! public routines
      public TestFlag,SetFlagParam,GetFlagTheta,GetFlagRadius

! public variables
    Logical, Public :: FlagActive
    
! Module Level Parameters

      
! Module Levels variables

    Real, Private, Save :: FlagRadiusMin
    Real, Private, Save :: FlagRadiusMax
    Real, Private, Save :: FlagTheta
contains

!*****************************************************************
! This is the routine use by part3 to test for the flag location

    Subroutine TestFlag(TH1,Theta,Y,E,nTurn,FSTP,RtnRegStp,ExitPart3)
! dummy variables
    real*8, intent(in) :: TH1 ! Angle of previous step
    real*8, intent(in) :: Theta ! Current Angle in radians
    real*8, intent(inout) :: Y(7)  !  - Y real*8 - vector of orbit values
    real*8, intent(in) :: E ! Energy
    integer, intent(in) :: nTurn ! Turn number
    real, intent (inout) :: FSTP ! Fractional Step
    Logical,intent(out) :: RtnRegStp !Signals a return to a regular step
    Logical,intent(out) :: ExitPart3 ! True to exit part3
! Local Variables
    Logical, Save :: Partial ! Executing the partial step to stripper
    Real*8, Save :: YLast(7) ! Store the orbit coordinates for partial step

    if (Partial) Then ! Means partial step is complete
        if(debug2) write(6,'("Flag Check TH,R",3f12.8)')theta*tcon,y(1),FlagRadiusMin,FlagRadiusMax
        if ((Y(1) .ge. FlagRadiusMin) .and. (Y(1) .le. FlagRadiusMax)) then 
            ExitPart3 = .true.
            write(6,'("***Stopped on flag***  Radius=",f12.3)')y(1)
            WRITE(55,'(" FLAG HIT AT",F10.3,6F9.3,11f9.3)')FlagTheta*tcon,Y(1),Y(2),Y(3)*tcon,E,Y(4),Y(5)
        else ! not yet at the foil radius - return to the previous location and continue integration
            Y = YLast
            RtnRegStp = .true.
        endif
        partial = .false.            
    else
        if ((TH1 .lt. FlagTheta) .and. (Theta .ge. FlagTheta))then ! flag angle was in the last step
            if(debug2) write(6,'("Flag possible,Th1,Th2,THflag ",6f12.8)')TH1*tcon,theta*tcon,FlagTheta*tcon
            partial = .true.
            Ylast = Y
            FSTP = (FlagTheta - TH1)/(Theta-TH1)
        else
            ExitPart3 = .false.
        endif
        Ylast = Y            
    endif

    end subroutine TestFlag
    
!************************************************
!  Set the Foil Parameters - called by parameter_load
!
    logical function SetFlagParam()

      real,external :: rdata_in
!      integer,external :: idata_in,cdata_in
      logical,external :: ldata_in    

        FlagActive = ldata_in(.false.,.false.)
        FlagTheta=rdata_in(.false.,0.0,360.,0.)/tcon
        FlagRadiusMin=rdata_in(.false.,0.0,1.E6,0.)
        FlagRadiusMax=rdata_in(.false.,FlagRadiusMin,1.E6,1.E6)
    
       SetFlagParam = .True.
    end function SetFlagParam
    
    real function GetFlagTheta()
        GetFlagTheta=FlagTheta*tcon
        return
    end function 
    real function GetFlagRadius()
        GetFlagRadius=FlagRadiusMin
        return
    end function
End module Flag