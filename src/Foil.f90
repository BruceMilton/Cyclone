module Foil
      use cyclone_lib
      use cyclone_data
      use bfields
      use io
      implicit none      
      private

! public routines
      public InitFoilRun, TestStripper,ManageStripper,WriteFoilDbs
      Public SetFoilParam

! public variables
    Logical, Public :: SuppressPrint
    Logical, Public :: FoilActive
    
! Module Level Parameters
    Integer,Private :: neq ! Max number of equations that can be integrated following stripper
    Parameter (neq = 17)
      
! Module Levels variables
    Logical, Private, Save :: Fine ! Select between course and fine search for stripper
    Logical, Private, Save :: EnergyDefined ! Selects between a radius or enegery defined foil
    Logical, Private, Save :: OptAngle ! True when stripper angle optimizing
    Logical, Private, Save :: TransferMatrices ! True to compute transfer matrix information
    Logical, Private, Save :: LogTransfer ! Output the transfer matrix info at each print step (unit 45)
! input values - should not change during a run
    Integer, Private, Save :: nFoilType=1 ! Select between simple,cross-over,energy
    Integer, Private, Save :: MaxIter=1

    Real, Private, Save :: Tolerence =1.e-5
    Real, Private, Save :: FoilRadiusInput
    Real, Private, Save :: FoilThetaInput
    Real, Private, Save :: FoilEnergy ! Energy to strip at when EnergyDefined = True
    Real, Private, Save :: FoilTilt = 0.0
    Real, Private, Save :: FoilDeltaE = 0.0 ! Energy loss in the foil
    Real, Private, Save :: StrippedCharge =1.00109 ! The charge after stripping - default is for H-
    Real, Private, Save :: CrossOverRadius
    Real, Private, Save :: CrossOverTheta
    Real, Private, Save :: StrippedStepSize = 2.0 !Step size after the stripper
    Integer, Private, Save :: nStrippedSteps = 100 ! number of steps after stripper
    Real, Private, Save :: FoilWidth ! radial range of foil
    Real, Private, Save :: FoilHeight ! vertical extent of foil
    Logical, Private, Save :: CheckFoilSize !if enabled then verify stripping position is in foil not frame
    Logical, Private, Save :: CheckFoilHeight !if enabled then verify the vertical stripping position is in a specific range
    Real, Private, Save :: FoilTop ! Upper valid value for z when CheckFoilHeight is true
    Real, Private, Save :: FoilBottom ! Lower valid value for z when CheckFoilHeight is true
    Real, Private, Save :: FoilExtent !Radial with of the foil if checking the height - defaults to very large

! Current active values
    Real, Private, Save :: FoilRadiusCurrent
    Real, Private, Save :: FoilThetaCurrent
    Real, Private, Save :: CourseTestRadius
    
! Values from the restart point
    Real*8, Private, Save :: ThetaStore,YStore(7),EStore
    Integer, Private, Save :: nTurnStore
        
contains
!*******************************************************
!  Initialize a Stripper run
!
    subroutine InitFoilRun
    
    if (debug2) write(6,'(" nFoilType = ",i3)')nFoilType
    select case (nFoilType)
    case (1) ! Fixed
        EnergyDefined = .false.
        OptAngle = .false.
    case (2) ! Cross - i.e. vary to hit defined cross-over
        EnergyDefined = .false.
        OptAngle = .true.
        ! prevent non-zero tilt
    case (3) ! E Fixed
        EnergyDefined = .true.
        OptAngle = .false.        
    case (4) ! E optimize
        EnergyDefined = .true.
        OptAngle = .true.
        ! prevent non-zero tilt
    case (5) ! Hold - use last run's foil position
        EnergyDefined = .false.
        OptAngle = .false.
        FoilRadiusInput = FoilRadiusCurrent
        FoilThetaInput = FoilThetaCurrent
        if(debug2) write(6,'("HOLD R,TH ",2f10.5)')FoilRadiusCurrent,FoilThetaCurrent
    case default
        write(6,'(" Illegal Foil Type")')
        stop
    end select
    
    FoilRadiusCurrent = FoilRadiusInput
    FoilThetaCurrent = FoilThetaInput
    CourseTestRadius = FoilRadiusInput ! There should be an input adjustment on this
    SuppressPrint = .False.
    
    if (optangle) Then ! setup storage
        !FoilThetaSave
        !CrossThetaSave
        !?PositionSave
    endif
    if (debug2) then
       write(6,'(" EnergyDefined,OptAngle ",2l5)')EnergyDefined,OptAngle
       write(6,'(" FoilRadiusCurrent,FoilThetaCurrent,CourseTestRadius ",3f12.5)')FoilRadiusCurrent,FoilThetaCurrent*tcon,CourseTestRadius
       write(6,'(" MaxIter,nStrippedSteps,Tolerence,StrippedStepSize ",2i5,2f10.4)')MaxIter,nStrippedSteps,Tolerence,StrippedStepSize
       write(6,'(" FoilEnergy,FoilDeltaE,StrippedCharge ",3f10.5)')FoilEnergy,FoilDeltaE,StrippedCharge
       write(6,'(" Transfer Matrices, TRN log ",2l5)')TransferMatrices,LogTransfer
       write(6,'(" Check foil vertical extent ",L5,2f10.3)')checkfoilheight,foiltop,foilbottom
    endif
    Fine = .false.
    return
    
    end subroutine InitFoilRun

!*********************
! Check for stripper location, and do integration if found

    Subroutine TestStripper(TH1,Theta,Y,E,nTurn,FSTP,RtnRegStp,ExitPart3)
! dummy variables
    real*8, intent(in) :: TH1 ! Angle of previous step
    real*8, intent(in) :: Theta ! Current Angle in radians
    real*8, intent(inout) :: Y(7)  !  - Y real*8 - vector of orbit values
    real*8, intent(in) :: E ! Energy
    integer, intent(in) :: nTurn ! Turn number
    real, intent (inout) :: FSTP ! Fractional Step
    Logical,intent(out) :: RtnRegStp !Singles a return to a regular step
    Logical,intent(out) :: ExitPart3 ! True to exit part3

! Local Variables
    Logical, Save :: Partial ! Executing the partial step to stripper
    Real*8, Save :: YLast(7) ! Store the orbit coordinates for partial step
    Real*8 :: deltaTheta ! change in theta from tilt
    Real*8 :: deltaR ! change if r for deltaTheta drift
    Real*8 :: PSQ
    Real :: FoilThetaTest ! Temp version of FoilThetaCurrent to be tested against last step location
    
    FSTP=0.0
    if(fine) then
!        if(debug2) write(6,'("Partial,Edef,TH1,Theta",2l5,2f10.3)')partial,energydefined,th1*tcon,theta*tcon
        if (Partial) Then ! Means partial step is complete
            if (EnergyDefined) Then
                FoilRadiusCurrent=Y(1)-1.E-5
                ExitPart3 = .true.
            else ! Radius
                if(debug2) write(6,'("Foil Check TH,R,Rfoil,Z",4f13.8)')theta*tcon,y(1),FoilRadiusCurrent,y(5)
                if (Y(1) .ge. FoilRadiusCurrent) then
                    if (CheckFoilHeight) then
                        if( y(4) .ge. FoilBottom .and. y(4) .le. FoilTop .and.Y(1) .lt. (FoilRadiusCurrent+FoilExtent)) then
                            ExitPart3 = .true.
                        else
                            Y = YLast
                            RtnRegStp = .true.
                            if(debug2) write(6,'("Failed vertical condition",5f13.8)')FoilBottom,y(4),FoilTop,y(1),FoilRadiusCurrent+FoilExtent
                        endif
                    else
                        ExitPart3 = .true.
                    endif
                else ! not yet at the foil radius - return to the previous location and continue integration
                    Y = YLast
                    RtnRegStp = .true.
                endif
            endif
            partial = .false.            
        else
            if( foilTilt .ne. 0)then ! calculate the foil tilt effect
              deltaTheta=foilTilt*(Y(1) - FoilRadiusCurrent)/Y(1) ! dt=tilt*(r-FoilCurrentRadius)
              FoilThetaTest=FoilThetaInput+deltaTheta
              deltaR=(FoilThetaTest-TH1)/(Theta-Th1)*(Y(1)-Ylast(1)) !estimate radius at foil angle
              deltaTheta=foilTilt*(Ylast(1)+deltaR - FoilRadiusCurrent)/(Ylast(1)+deltaR) ! re-calculate the foil angle at new radius
              FoilThetaTest=FoilThetaInput+deltaTheta
            else
              FoilThetaTest=FoilThetaCurrent
            endif
            if ((TH1 .lt. FoilThetaTest) .and. (Theta .ge. FoilThetaTest))then ! foil was in the last step
                if(debug2) write(6,'("Foil possible,Th1,Th2,THfoil ",6f13.8)')TH1*tcon,theta*tcon,FoilThetaTest*tcon
                if(debug2) write(6,'("r,rlast,deltaR,Rest",6f13.8)')y(1),ylast(1),deltaR,ylast(1)+deltaR
                partial = .true.
                Ylast = Y
                FoilThetaCurrent=FoilThetaTest
                FSTP = (FoilThetaCurrent - TH1)/(Theta-TH1)
            else
                ExitPart3 = .false.
            endif
            Ylast = Y            
        endif
    else ! Still in course search phase
        if (EnergyDefined) then
            if(E .ge. FoilEnergy) then
                call StartFine(Theta,Y,E,nTurn)
                ExitPart3 = .true.
                return
            else
                ExitPart3 = .false.
                return
            endif
        else
            if (y(1) .gt. CourseTestRadius) then
                call StartFine(Theta,Y,E,nTurn)
                ExitPart3 = .true.
            else
                ExitPart3 = .false.
                return
            endif
        endif        
    endif

    end subroutine TestStripper

!*********************************************************
! This routine controls the various intergration steps that follow a course check identification of the
!   stripper    
!
    subroutine ManageStripper(xonly)

!*******************************************************    

! Dummy Variables
    Logical, intent(in) :: xonly
! Local Variables
    Integer :: iter
    Logical,  Save :: Done ! 
    Logical,  Save :: Short ! Did not make cross-over on last run
    Logical,  Save :: FinalIter ! Final iteration
    Real*8 :: YI(7),EI
    Real :: THI,THF
    Integer :: nTurn,ith_sw
    Real :: ThetaCross ! This the interpolated THETA for the Cross-Over Radius (CrossOverRadius)
    Real :: FoilThetaPrevious,ThetaCrossPrevious ! These maintain the previous run's locations
    Real :: Tmag,Err ! used when optimizting the foil location
    Real :: StrippedRadius
    
    iter = 1
    Done = .false.
    if (OptAngle) then
        FinalIter = .false.
        SuppressPrint = .true.
        write(6,'(" #ITER   INITIAL R  INITIAL TH   FINAL R    FINAL TH     TH MAG.   TH ERR(DEG)")')
    else
        FinalIter = .true.
        SuppressPrint = .false.
    endif
    if(debug3) SuppressPrint=.false.
    ThetaCrossPrevious=0.0
    FoilThetaPrevious=0.0

    do While (.Not. Done)
        YI=YStore
        EI=EStore
        THI = ThetaStore*tcon
        nTurn = nTurnStore
        ith_sw = 0
        call dummyPart3(YI,EI,nTurn,THI,ith_sw,1.0,xonly,*999,*500)
        if (ExitCode .ne. iex_FlagHit) then
            !stop "Exited Part 3 not stripped"
            ExitCode = iex_nostrip
            write(6,'("Exiting Part 3 not stripped")')
        endif
        return
500     continue
        if (FinalIter) then
            if (CheckFoilSize) then
              if((YI(1) .gt. (FoilRadiusCurrent + FoilWidth)) .or. abs(YI(4)) .gt. FoilHeight) then
                Write(6,'(" ")')
                write(6,'("Particle hit foil frame,r,z=",2f10.2,"(",2f10.2,")")')YI(1),YI(4),FoilRadiusCurrent + FoilWidth,FoilHeight
                write(59,'(F10.3,6F9.3,11f9.3)')THI,YI(1),YI(2),YI(3)*tcon,EI,YI(4),YI(5)
                ExitCode=iex_FoilFrame
                return
              endif
            endif
            Write(6,'(" ")')
            write(6,'("******  STRIPPING OCCURS AT R= ",F10.5,"  THETA= ",F10.5," ********")')YI(1),THI
            ! if tilt is non-zero then calc error in angle and print
            if(FoilTilt .ne. 0) then
              THF=(FoilThetaInput+FoilTilt*(YI(1)-FoilRadiusCurrent)/YI(1))*TCON
              write(6,'("Foil Tilt correction, actual=",f10.5," Diff=",1PE12.4)')THF,(THF-THI)
            endif
            write(6,'(" Nturn  Theta      R       Pr      PTheta    E        S")')
            WRITE(56,'(" STRIPPING OCCURS AT",F10.3,6F9.3,I7)')THI,YI(1),YI(2),YI(3)*tcon,EI,YI(4),YI(5),nrun
                        ExitCode = iex_Stripped
        endif
        ! Capture the actual stripping radius (always > FoilRadius)
        StrippedRadius = YI(1)
        ! Update the foil angle with the angle found
        FoilThetaCurrent = THI/tcon
        ! Now calculate the trajectory from the foil outwards
        if (.Not. StripTraj(YI,EI,nTurn,THI,ThetaCross,FinalIter,Done,xonly))then ! failure condition
            write(6,'(" Strip Trajectory Failed ")')
            return
        endif
        if(OptAngle .and. (.not. FinalIter)) then
            Tmag = 1.0
            Err = CrossOverTheta - ThetaCross
            if (Iter .gt. 1) then
                Tmag = (ThetaCross - ThetaCrossPrevious)/(FoilThetaCurrent - FoilThetaPrevious)
                if (abs(tmag) .lt. 0.25) TMAG = Sign(0.25,Tmag)
            endif
            Write(6,'(I5,8f12.5)')ITER, StrippedRadius,FoilThetaCurrent*tcon,CrossOverRadius,ThetaCross*tcon,Tmag,err*tcon
            ! Update the stored locations for the next run
            ThetaCrossPrevious = ThetaCross
            FoilThetaPrevious = FoilThetaCurrent
            ! Calculate the new stripper location
            FoilThetaCurrent = FoilThetaCurrent + Err/Tmag
            if (abs(Err) .lt. tolerence) Then
                FinalIter = .true.
            endif
        endif
        iter = iter +1
        if (iter >= MaxIter) FinalIter = .true.
    end do
    
    if(logtransfer) write(45,*) ! leave a blank line between orbits
    
    return
999 Write(6,'(" Part 3 in stripper error return")')
    ExitCode = iex_nostrip
    return
    end subroutine ManageStripper

!*******************************************************
! Store values for restarting the particle
!
    Subroutine StartFine(Theta,Y,E,nTurn)
!*******************************************************
! dummy variables
    real*8, intent(in) :: Theta ! Current Angle in radians
    real*8, intent(in) :: Y(7)  !  - Y real*8 - vector of orbit values
    real*8, intent(in) :: E ! Energy
    integer,intent(in) :: nTurn

! local variables
    real :: tau
    
    ThetaStore = Theta
    YStore = Y
    EStore = E
    nTurnStore = nTurn
    Fine = .true.
    if (OptAngle) SuppressPrint = .true. ! turn off Part3 printing
    if (debug2) write(6,'(" Initializing a fine search, Theta = ",2f10.5)')theta,YStore(6)
    ! write restart data to a file
    TAU=Y(3)
    CALL THNRM(tau)
    TAU=TAU*tcon    
    write(58,'("THETA=",f9.3,",ANGLE, ENERGY=",E14.7," R=",E14.7," PR=",E14.7," TAU=",E14.7,/," Z=",E14.7," PZ=",E14.7," RUN")')Theta*tcon,E,Y(1:2),tau,Y(4:5)
        
    end subroutine StartFine
        
!**********************************************************************************
! Integrate the particle trajectory after the stripper
!
    Logical Function StripTraj(Y,E,nTurn,Theta,ThetaCross,FinalIter,Done,Xonly)
! dummy variables
      real*8, intent(inout) :: Y(7)  !  - Y real*8 - vector of orbit values
      Real*8, intent(in) :: E
      Integer, intent(in) :: nTurn
      Real, intent(in) :: Theta
      real, intent(out) :: ThetaCross ! This is the interpolated THETA for the Cross-Over Radius (CrossOverRadius)
      logical, intent(inout) :: FinalIter ! Indicates that we are now ready for last iteration
      logical, intent(inout) :: done ! Indicated have completed last interation
      logical,intent(in) :: xonly ! true is on radial motion, false means vertical motion also
      
! local variables
      integer :: neq_used ! Actual active number of equations of motion

      real*8 RRK,CK(neq),Q(neq),YS(neq)
      real*8 ARK(4)/.5D+0,.292893219D+0,1.707106781D+0,.1666666667D+0/
      real*8 BRK(4)/2.D+0,1.D+0,1.D+0,2.D+0/
      real*8 CRK(4)/-.5D+0,-.292893219D+0,-1.707106781D+0,-.5D+0/
      
      Real*8 :: PSQ,GAM,P,ES
      Real*8 :: temp1,temp2,temp3
      Real*8 :: Step ! The integration step size

      logical :: Loop
      logical :: CrossOverFound
      Integer :: istp,i,j
! Save the last three R,Th values
      Real*8 :: CoordStore(neq,2)
      real :: Outputs(neq)
      
! Convience mapping of intergartion variables
      real*8 ::R,PR,PTHETA,S,Z,PZ
      EQUIVALENCE (R,YS(2)),(PR,YS(3)),(PTHETA,YS(4)),(S,YS(5)),(Z,YS(6)),(PZ,YS(7))
      real*8 :: Z1,PZ1,Z2,PZ2
      EQUIVALENCE (Z1,YS(8)),(PZ1,YS(9)),(Z2,YS(10)),(PZ2,YS(11))
      real*8 :: X1,PX1,X2,PX2
      EQUIVALENCE (X1,YS(12)),(PX1,YS(13)),(X2,YS(14)),(PX2,YS(15))
      real*8 :: D1,D2
      EQUIVALENCE (D1,YS(16)),(D2,YS(17))

! Figure out if transfer matrices are used or not
        if (TransferMatrices) then
            neq_used = neq
        else
            neq_used = 7
        endif
! Setup up the integration variables
        Step = StrippedStepSize
        Q=0
        YS(1) = Theta/tcon
        YS(2) = y(1) ! R
        YS(3) = y(2) ! PR
        ES = E - FoilDeltaE
        PSQ=ES/E0
        GAM=1.D+0+PSQ
        PSQ=PSQ*(2.D+0+PSQ)*ASQ
        P = sqrt(psq)
        YS(4) = Sqrt(PSQ - y(2)**2) ! Ptheta
        YS(5) = 0.0 ! S
        YS(6) = Y(4) !z1
        YS(7) = Y(5) ! PZ1
        ! Note Y(6) is stripping loses
        ! Transfer matrix equations
        YS(8) = 1.0 !z1
        YS(9) = 0.0 ! PZ1
        YS(10) = 0.0 !Z2
        YS(11) = P/1000 ! PZ2
        YS(12) = 1.0 ! x1
        YS(13) = 0.0 ! x1' 
        YS(14) = 0.0 ! x2
        YS(15) = 0.001 !x2'
        YS(16) = 0.0 !d1
        YS(17) = 0.0 !d2
        outputs = ys
        outputs(9)=pz1*1000/p
        outputs(11)=pz2*1000/p
        outputs(13)=px1*1000
        outputs(15)=px2*1000
        outputs(16)=(x2*D1-x1*D2)*1000.
        outputs(17)=(px2*D1-px1*D2)*1.D6                

        if (FinalIter) Then
!            write(6,200)nturn,ys(1)*tcon,' ',(ys(j),j=2,4),ES,YS(5),YS(6),YS(7)
            if (LogTransfer) then
                write(45,201)ys(1)*tcon,ys(2),(outputs(j),j=8,17)
            endif
            call DoPrint(ys,es,nturn,xonly)
        endif
        
! Initialize the storage to the starting position
      do i = 1,2
        CoordStore(1:neq,i) = outputs(1:neq)
      enddo
! Init Switches
      CrossOverFound = .False.
      loop = .true.
      istp = 0

      do while (loop)
        istp=istp+1
        do j = 1,4
            IF (.not. BFLDCAL(YS(1),YS(2),.false.,1) ) go to 113
    !       Compute orbit derivatives
            CK(1)=STEP*YS(4)/(YS(2)*P) ! d theta
            CK(2)=STEP*YS(3)/P         ! d r
            CK(3)=CK(1)*(YS(4)+StrippedCharge*YS(2)*bz) !d p_r
            CK(4)=CK(2)*(-YS(4)/YS(2)-StrippedCharge*bz)    !d p_theta
            CK(5)=STEP !ds
    !       vertical equations for ray
            CK(6)=STEP*pz/P 
            CK(7)=-CK(1)*StrippedCharge*R*Z*(dBZdR-Pr*dBZdT/(R*Ptheta))
    !       vertical equations dz1,dpz1,dz2,dpz2 of matrix
            CK(8)=STEP*pz1/P 
            CK(9)=-CK(1)*StrippedCharge*R*Z1*(dBZdR-Pr*dBZdT/(R*Ptheta))
            CK(10)=STEP*pz2/P
            CK(11)=-CK(1)*StrippedCharge*R*Z2*(dBZdR-Pr*dBZdT/(R*Ptheta))
    !       transfer matrix equations
            temp1=step*StrippedCharge*(dbzdr*YS(4)-YS(3)*dBZdT/YS(2))/psq
            temp2=-step*(StrippedCharge*bz)**2/psq
            temp3=-step*(StrippedCharge*bz)/p
            ck(12)=step*YS(13) ! x1
            ck(14)=step*YS(15) ! X2
            ck(13)=YS(12)*(temp1+temp2) !x1'
            ck(15)=YS(14)*(temp1+temp2) !x2'
            ck(16)=temp3*YS(12) !D dispersion
            ck(17)=temp3*YS(14) !D' dispersion
            do i=1,neq_used
              RRK=ARK(J)*(CK(I)-BRK(J)*Q(I))
              YS(I)=YS(I)+RRK
              Q(I)=Q(I)+3.D+0*RRK+CRK(J)*CK(I)
            end do  
        end do
        outputs = ys
        outputs(9)=pz1*1000/p
        outputs(11)=pz2*1000/p
        outputs(13)=px1*1000
        outputs(15)=px2*1000
        outputs(16)=(x2*D1-x1*D2)*1000.
        outputs(17)=(px2*D1-px1*D2)*1.D6                
        CoordStore(1:neq_used,1)=CoordStore(1:neq_used,2)
!        CoordStore(1:neq_used,2)= YS(1:neq_used)
        CoordStore(1:neq_used,2)= outputs(1:neq_used)
        if (.Not. CrossOverFound .and. (ys(2) .ge. CrossOverRadius)) then ! past x-over pt
          ! calc theta
          CrossOverFound = .true.
          Call FindTheta(CoordStore,CrossOverRadius,ThetaCross)
          if (FinalIter  .or. Debug3) then
            call FindAtTheta(CoordStore,ThetaCross,Outputs,neq_used)
            Outputs(1)=Outputs(1)*tcon
            if (xonly) then
                write(6,200)nturn,Outputs(1),'CO',(Outputs(j),j=2,4),ES,Outputs(5)
            else
                write(6,200)nturn,Outputs(1),'CO',(Outputs(j),j=2,4),ES,(Outputs(j),j=5,7)
            endif
            write(57,'(" Cross Over Point At",F10.3,6F9.3,11f9.3)')(Outputs(j),j=1,4),ES,(Outputs(j),j=5,neq_used)
          endif
          if (.Not. FinalIter) loop = .false.
        endif
        if (FinalIter .or. Debug3) then
            call DoPrint(ys,es,nturn,xonly)
            if (LogTransfer) then
                write(45,201)ys(1)*tcon,ys(2),(outputs(j),j=8,17)
            endif            
        endif
        if (istp .ge. nStrippedSteps) loop = .false.
      end do
      
      if( .Not. CrossOverFound) then
            write(6,'(" **** Failed to find Cross-Over ****")')
      endif
      if (finalIter) done = .true.
      StripTraj = .true.
      return
      
113   if (CrossOverFound) Then
        if (finalIter) done = .true.
        StripTraj = .true.
      else
        Write(6,'("Magnetic Field Failure ")')
        StripTraj = .false.
      endif
      return
  200   FORMAT(' ',I4,F8.2,A2,11F9.4)        
  201   Format(" ",F8.2,1X,11F9.4)
    end function

!************************************************
!  Manage output from StripTraj
    Subroutine DoPrint(ys,es,nturn,xonly)
!************************************************
! dummy variables
      real*8, intent(in) :: YS(*)  !  - Y real*8 - vector of orbit values
      Real*8, intent(in) :: ES
      Integer, intent(in) :: nTurn
      logical,intent(in) :: xonly ! true is on radial motion, false means vertical motion also      
 
 ! local variables
      integer :: j
           
        if(xonly) then
            write(6,200)nturn,ys(1)*tcon,' ',(ys(j),j=2,4),ES,YS(5)
        Else
            write(6,200)nturn,ys(1)*tcon,' ',(ys(j),j=2,4),ES,(YS(j),j=5,7)
        endif
        If (io_control35.active) then
            io_control35.lines=io_control35.lines+1
            if (io_control35.formatted) then
                WRITE(35,'(9e14.6,i2)')ys(1)*tcon,ys(2),ys(3),0.0,ES,0.0,ys(6),ys(7),nturn,4
            else
                write(35)sngl(ys(1)*tcon),sngl(ys(2)),sngl(ys(3)),0.0,sngl(es),0.0,sngl(ys(6)),sngl(ys(7)),sngl(nturn),4.
            endif
        endif
         return
 
  200   FORMAT(' ',I4,F8.2,A2,11F9.4)        
  201   Format(" ",F8.2,1X,11F9.4)
    end subroutine

!************************************************
!  Interpolate to find Theta for a specified Radius
    Subroutine FindTheta(Coord,Rin,Tout)
!************************************************
        Real*8, intent(in) :: Coord(neq,2)
        Real, intent(in) :: Rin
        Real, intent(out) :: Tout

        Tout = Coord(1,2) + (Coord(1,1) - Coord(1,2))*(Rin-Coord(2,2))/(Coord(2,1)-Coord(2,2))
    end subroutine FindTheta

!************************************************
!  Interpolate to find values at a specified Theta Value
    Subroutine FindAtTheta(Coord,Theta,Output,n)
!************************************************
        Real*8, intent(in) :: Coord(neq,2) ! Contains the orbit values before and after interpolation point
        Real, intent(in) :: Theta
        Integer, intent(in) :: n ! number for active equations

        Integer ::i
        Real :: FST
        Real :: Output(neq)
        
        FST = (theta-Coord(1,1))/(Coord(1,2) - Coord(1,1))
        do i = 1 ,n
            output(i) = Coord(i,1) + fst *(Coord(i,2)-Coord(i,1))
        enddo

    end subroutine FindAtTheta

!************************************************
!  Set the Foil Parameters
!
    logical function SetFoilParam()
      character*10 foiltypes(5)/'FIXED','CROSS','EFIXED','EOPT','HOLD'/ ! Should change 'cross' as a type!!!
      character*10 FoilCmds(8)/'OPT','CROSS','ACTIVE','POSITION','ENERGY','TRANSFER','SIZE','HEIGHT'/
      integer iFoilCmd
      
      real,external :: rdata_in
      integer,external :: idata_in,cdata_in
      logical,external :: ldata_in

        iFoilCmd = cdata_in(.false.,FoilCmds,8,1)
        select case (iFoilCmd)
        case (1) ! OPT
            nFoilType = cdata_in(.false.,Foiltypes,5,1)
            MaxIter = idata_in(.false.,1,100,10)
            Tolerence = rdata_in(.false.,1.e-8,1.0,0.001)/tcon
            if(debug2) write(6,'("OPT ",3I5)')nFoilType,MaxIter,Tolerence
        case (2) ! Cross
            CrossOverRadius = rdata_in(.false.,0.,10000.,0.)
            CrossOverTheta = rdata_in(.false.,0.,360.,0.)/tcon
        case (3) ! Active
            FoilActive = ldata_in(.false.,.true.)
            nStrippedSteps = idata_in(.false.,1,10000,100)
            StrippedStepSize=rdata_in(.false.,0.0,1000.,2.)
            if(FoilActive .and. FoilRadiusInput .le. 0.0) then
                FoilActive = .false.
                write(6,'(" Must Specify foil position before activating")')
            endif
            if(debug2) write(6,'("Active ",l3,i5,f10.3)')FoilActive,nStrippedSteps,StrippedStepSize
        case (4) ! Position
            FoilRadiusInput =rdata_in(.false.,0.,10000.,0.)
            FoilThetaInput = rdata_in(.false.,0.,360.,0.)/tcon
            FoilTilt = rdata_in(.false.,-90.,90.,0.)/tcon
            FoilActive = .true.
        case (5) ! Energy
            FoilEnergy = rdata_in(.false.,0.,10000.0,0.)
            FoilDeltaE = rdata_in(.false.,-1.,1.0,0.)
            StrippedCharge = rdata_in(.false.,-100.,100.,1.00109)
        Case (6) !Transfer
            TransferMatrices = ldata_in(.false.,.true.)
            LogTransfer = ldata_in(.false.,.false.)
            if (.not. TransferMatrices) LogTransfer = .false. ! can't log if not calculated
            if (LogTransfer) then
              if (openFile(45,return_msg) .ne. 0) then
                SetFoilParam=.false.
                return
              endif
            endif
        case (7) ! Size = Frame Check
            CheckFoilSize = ldata_in(.false.,.true.)
            FoilWidth = rdata_in(.false.,0.,1000.0,1.)
            FoilHeight = rdata_in(.false.,0.,1000.0,1.)
        case (8) ! Height = Foil has limted vertical extent
            CheckFoilHeight = ldata_in(.false.,.true.)
            FoilTop = rdata_in(.false.,-1000.,1000.,1.)
            FoilBottom = rdata_in(.false.,-1000.,1000.,0.)
            if (FoilTop .lt. FoilBottom) write(6,'("****Warning the foil top is below the foil bottom so zero extent *****")')
            FoilExtent = rdata_in(.false.,-1000.,1000.,1000.)
        case default
            setfoilparam = .false.
            Write(6,'('' Illegal foil command '')')
            return
        end select
! need to add a radius offset for search start
! Also need to account for back bend particles during search
        SetFoilParam = .true.
    end function SetFoilParam

!****************************************************************
! write foil information on dbs (unit 34 file)
!    
    subroutine WriteFoilDbs(iunit,iline)
      integer, intent(in) :: iunit ! unit to write dbs data
      integer, intent(in) :: iline ! card number - normally 6
      
      real,external :: rout
      
      write(iunit,326)iline,rout(foilActive),float(nFoilType),CrossOverRadius,CrossOverTheta*tcon,FoilRadiusInput, &
        FoilThetaInput*tcon,FoilTilt*tcon,FoilEnergy
 326  format(i3,8f9.4)
      
    end subroutine
    

end module foil