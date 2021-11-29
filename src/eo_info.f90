! Contains the equilibrium orbit routines
! The public routines do the following
!   eo_input - read EO data from a file
!   eo_output - write EO data to file (for future read)
!   eo_make - calculate EO, info is written to EO.log
!   load_ellipse - return orbit coordinates from a stored eigen-ellipse
!

! sigma= SIGMA is compiler variable used to link in the libraries, 0 to debug without lib, 1 is with
    
      include 'mkl_vsl.f90'
      module eo_info

      use cyclone_data
      use cyclone_lib
      use bfields
      use IO
#if SIGMA == 1
      USE MKL_VSL_TYPE, ONLY:VSL_STREAM_STATE, VSL_BRNG_MCG31, VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2, VSL_MATRIX_STORAGE_FULL
      USE MKL_VSL, ONLY:vslnewstream, vsrnggaussianmv ! spotrf
#endif

      implicit none
      private

! public routines in this module
      public eo_input,eo_output,eo_make,load_ellipse_info,fill_ellipse
      public SetEllipseFillType, SetEllipseNumber, SetEllipseOffset, SetEllipseType
      public SetVaryTau,GetPot,Get_EO,SetEllipseBound, ReadEllipseOffset
      public clean_EO_info

!public variables

      INTEGER,PARAMETER,public :: NEO_MAX=1000   !Maximum number of eo values
      real,public,save,allocatable :: reop(:,:),preop(:,:) ! storage for all theta EOs
      real,public,save,allocatable :: REO(:,:,:) ! first index 1=r, 2=pr, 2nd index=ispr
      real,public,save :: EOI=0.0,DEO=0.0,THEO=-1.0
      integer,public,save :: NEO=0,ieo_type
      integer,public,save,allocatable :: EO_NUMBER(:) 
      integer,public,save :: ieo_stp=-2 ! step number for fixed theo eos
      integer,public,save :: istp_nsc ! number of steps/sector for fixed theta
      integer,public,save :: kstep_eo ! numer of EOs/energy - all theta eo
      integer,public,save :: eo_ngap ! number of gap locations
      real,public,save :: ellipse_area(2) ! eigen ellipse info
      integer,public,save :: ellipse_type ! eigen ellipse info 1=x,2=z
      integer,public,save :: ellipse_counter ! counter for ellipse particles
      logical,public,save :: run_ellipse ! particles start on ellipse
      logical,public,save :: ellipse_norm ! true if area is normalized
      
! Module level variables
      integer,private,save :: nEllipseParticles ! Number of particles the will be used to fill/outline ellipse
      integer,private,save :: nEllipseFilltype ! Type of Ellipse fill - defined in parameter_load
      real*8,private,save :: Ellipse_Offset_R ! the ellipse centroid r value or XI
      real*8,private,save :: Ellipse_Offset_PR ! the ellipse centroid pr value or PXI
      real,private,save :: Ellipse_Offset_Z ! the ellipse centroid z value or ETA
      real,private,save :: Ellipse_Offset_PZ ! the ellipse centroid pz value or PETA
      logical,private,save :: UseEllipseOffSet ! False means use the EO values
      real(kind=4),private,save,allocatable :: ellipse(:,:) ! Particle Coordinates for ellipses
      real(kind=4),private,save,allocatable :: tau(:) !Tau distribution for ellipses if enabled
      real,private,save :: sigma_in(10) ! if sigma values entered, save here
      logical,private,save :: UseEigen ! Use the eigen ellipse values for sigma
      logical,private,save :: UseCorr=.false. ! Indicates that x-z correlations input
      logical,private,save :: VaryTau=.false. ! Switch to enable randomizing Tau values
      real,private,save :: TauMin !Lower limit of Tau values
      real,private,save :: TauRange ! delta Tau
      real,private,save:: VaryPot ! this determines if the VaryTau is Tau or Phi
      real,private,save:: AreaScale ! this converts from mm-mrad to length-length
      real,private,save:: pScale ! this converts from mrad to length
      real*8,private,save:: y(12) !EO integration variables
      real,private,save:: arg(2) !the ellipse angle argument, orginally dimensioned in make_eo
      real,private,save:: GaussianBound ! trancation boundary of gaussian dist in terms of simgas
      logical,private,save:: GaussTruncate=.false. ! if trancuation is enable
#if SIGMA == 1
      TYPE (VSL_STREAM_STATE) :: stream
#endif
      Logical,private,save :: StreamInitialized=.false.
      
      contains

!***********************************************************

      subroutine eo_input(name)
      CHARACTER*csize NAME,name44
      integer :: ichar,i,j,k,ierr
      integer :: nlist ! input variable from unit 25 - should be 4*ndee

      OPEN (UNIT=25,STATUS='OLD',READONLY,file=name,iostat=ierr)
      if(ierr .ne. 0) then
         write(6,'(''***Error Opening EO file ****'',i5)')ierr
         write(6,'(a<csize>)')name
         neo=0
         return
      endif
      inquire(25,name=name44)
      name44=adjustl(name44) !remove leading blanks
      ichar=len_trim(name44)
      write(6,'('' EO Data from: '',a<ichar>)')name44(1:ichar)

! first line gives file details

      read(25,'(2i10,3e16.8)',iostat=ierr)ieo_type,neo,eoi,deo
      if(ierr .ne. 0) then
         write(6,'(''***Error Reading EO file header ****'',i5)')ierr
         return
      endif

! ieo_type selects type of eo data
!          =1 not supported
!          =2 at one fixed angle th=theo - part2 
!          =3 at gap crossings - part3 - LIKE SPRGAPZ
!          =4 at all angles

! EO input
!


      IF(NEO .GT. neo_max) STOP 'TOO MUCH EO DATA'
      IF(NEO .LT. 1) then
         write(6,'(''*** NEO=0 ***'')')
         neo=0
         return
      endif
      IF(DEO .LE. 0.0) then
         write(6,'(''*** DEO=0 ***'')')
         neo=0
         return
      endif

      if(ieo_type .le. 2) then
         read(25,'(2i5)',iostat=ierr)ieo_stp,istp_nsc
         if(allocated(reo)) deallocate(reo)
         allocate(reo(2,1,neo),stat=ierr)
         call allocate_check(ierr,'EO Storage reo')
         READ(25,'(2f12.5)',iostat=ierr)(REO(1,1,I),REO(2,1,I),I=1,NEO)
         if(ierr .ne. 0) then
            write(6,'(''***Error Reading EO file ****'',i5)')ierr
            stop 'EO File Error'
         endif
      else if(ieo_type .eq. 3) then
        if(ndee .le. 0) then
           write(6,'('' Must issue DEE command before reading type 3 EO data '')')
           stop 'bad Input'
        endif
        read(25,'(2i5)',iostat=ierr)nlist,eo_ngap
        if(nlist .ne. 4*ndee) then
           write(6,'('' Warning Incorrect Gap Specs on EO File '')')
           nlist=min(nlist,4*ndee)
        endif
        if(eo_ngap .lt. 1) then
           write(6,'(''*** Must be at least 1 gap in EO file'')')
           stop 'EO File Error'
        endif
        if(allocated(EO_NUMBER).and.(size(EO_NUMBER).ne.4*ndee)) deallocate(EO_NUMBER)
        if(.not. allocated(EO_NUMBER)) allocate (EO_NUMBER(4*ndee),stat=ierr)
        call allocate_check(ierr,'Gap Array EO_NUMBER')
        read(25,'(20i5)')(EO_NUMBER(i),i=1,nlist)
        istp_nsc=1
        if(allocated(reo)) deallocate(reo)
        allocate(reo(2,eo_ngap,neo),stat=ierr)
        call allocate_check(ierr,'EO Storage reo')
        do i=1,neo
           read(25,'(8F12.5)',iostat=ierr)((reo(k,j,i),k=1,2),j=1,4)
           if(ierr .ne. 0) then
              write(6,'(''***Error Reading EO file ****'',i5)')ierr
              write(6,'('' NEO = '')')i
              stop 'EO File Error'
           endif
        end do
      elseif(ieo_type .eq. 4) then
        read(25,'(2i5)',iostat=ierr)kstep_eo,istp_nsc
        do i=1,neo
            read(25,'(8F12.5)',iostat=ierr)(reop(j,i),preop(j,i),j=1&
                 & ,kstep_eo)
           if(ierr .ne. 0) then
              write(6,'(''***Error Reading EO file ****'',i5)')ierr
              write(6,'('' NEO = '')')i
              stop 'EO File Error'
           endif
        end do
      else
         write(6,'(''*** Illegal EO type ***'',i5)')ieo_type
         stop 'EO File ERROR'
      endif
      return
      END subroutine eo_input

!***********************************************************
      subroutine eo_output(name)

      integer :: I,j,l,ierr

      CHARACTER*csize NAME

      !OPEN (UNIT=26,STATUS='UNKNOWN',file=name,iostat=ierr)
      ierr=openFile(26,return_msg)
      if(ierr .ne. 0) then
         write(6,'(''***Error Opening EO Output file ****'',i5)')ierr
         write(6,'(a<csize>)')name
         return
      endif

      IF(NEO .lt. 1) Then
         write(6,'('' No EO data to write'')')
      endif

      write(26,'(2i10,1p,3e16.8)')ieo_type,neo,eoi,deo,theo

      if(ieo_type .le. 2) then
         write(26,'(2i5)')ieo_stp,istp_nsc
         DO I=1,NEO
            WRITE(26,'(2f12.5)')(REO(1,1,I),REO(2,1,I))
        enddo
      else if(ieo_type .eq. 3) then
        write(26,'(2i5)')ndee*4,eo_ngap
        write(26,'(20i5)')(EO_NUMBER(i),i=1,ndee*4)
        do i=1,neo
           write(26,'(8F12.5)')((reo(l,j,i),l=1,2),j=1,eo_ngap)
        enddo
      elseif(ieo_type .eq. 4) then
        write(26,'(2i5)')kstep_eo,istp_nsc
        do i=1,neo
           write(26,'(8F12.5)')(reop(j,i),preop(j,i),j=1,kstep_eo)
        end do
      endif
      close(26)
      return
      END subroutine eo_output
!***********************************************************
!  Calculate EO at the defined energies (neo,deo,eoi)
!   EO data is output to EO.log (unit 20), and this include nur, nuz & w0/w-1
!   EO data is also collected at specific angles depending on the value or ieo_type
!       1 = Used to calculate EO & possibly its eigen-ellipse (that is stored & then recouvered using load_ellipse)
!       2 = At a fixed angle, defined as THEO (for part2)
!       3 = At the gaps (same as SPRGAPZ)
!       4 = At all print angles
!************************************************************
      subroutine eo_make(rguess,prguess)

        use spiral_gaps

! dummy variables
        real*4 :: rguess,prguess

! local variables
      integer :: I,j,k,i1,ierr
      real*8 :: peo

      real*8 :: PSQ,gam,q(12),pts,pt
      real*8 :: rg,prg ! guesses for the start of the next iteration
      real*8,allocatable ::rst(:),prst(:) ! temp storage arrays
      integer :: IEO ! main energy loop counter
      integer :: nstp ! number of steps for integration
      integer :: istp ! loop counter for integration
      real*8 :: STP !step size in radians
      real*8 :: SECT ! radians per sector = integration angle
      integer :: ith ! magnetic field position
      logical :: LastIter ! set true for last iteration of eo search
      integer :: ntry ! search loop counter
      Integer :: nxcros,nzcros ! crossing counters
      real  :: yp5,yp10 ! for tracking crossings
      real  :: sn,cn ! sin and cos terms for focusing freq
      real  :: rznt,rznu(2) ! nr and nz
      real*8 :: RRK,CK(12) !RRK variables ck=derivatives
      REAL*8  ARK(4)/.5D+0,.292893219D+0,1.707106781D+0,.1666666667D+0/
      real*8  BRK(4)/2.D+0,1.D+0,1.D+0,2.D+0/
      real*8  CRK(4)/-.5D+0,-.292893219D+0,-1.707106781D+0,-.5D+0/
      integer :: neq ! number of integration variables
      real*8 :: EP1,EP2,DEN !convergance variables
      real*8 :: xc0,xc1,xc2,xc3,xc4 ! intermediate results
      real*8 :: PHS,PHSP,FOE ! phase variables
      real*8 :: RAV ! average radius
      real   :: E ! energy of eo
      logical :: bad ! set true if an EO is not found
      integer :: igp ! gap loop counter
      real :: fth ! fraction of step
      real :: dummy ! for call to thin = spiral constant
      logical :: found ! indicates a gap has been found
      real,allocatable :: cross(:) ! temp storage of crossing angles


      NSTP=NTHETA/2
      SECT=TPI/NSC
      STP=SECT/NSTP
      istp_nsc=nstp

      write(6,'('' Computing EOs n,ei,de'',i5,2f10.5)')neo,eoi,deo
      if(ieo_type .ne. 1) then
         !OPEN(26,status='unknown',iostat=ierr)
         ierr=openFile(26,return_msg)
         if(ierr .ne. 0) then
            write(6,'(''Error opening unit 26 for EOs '',i5)')ierr
            stop 'file error'
         endif

        write(26,308)E0,chg,nsc,acon,bcon,anu0,nint(sect*tcon),ndee
 308    FORMAT(5H E0 =,F10.3,5X,5HQ/E =,F9.6,I6,8H SECTORS// &
             4H A =,F10.5,5X,3HB =,F9.5,5X,5HNU0 =,F9.5// &
             22H INTEGRATION FROM 0 TO,I4,8H DEGREES,I6,5H DEES/)
        write(26,316)
 316    FORMAT(4X,1HE,9X,3HRAV,7X,4HR(0),6X,5HPR(0),4X,6HW0/W-1,5X, &
             4HNU R,6X,4HNU Z,8X,4HF(E))
      endif

! allocate temp storage

      allocate(rst(nstp+1),stat=ierr)
      call allocate_check(ierr,'EO Storage - rst')
      allocate(prst(nstp+1),stat=ierr)
      call allocate_check(ierr,'EO Storage -prst ')


      RG=rguess
      PRG=prguess
      if(RG .le. 0.0) then
         GAM=EOI/E0+1.D0
         PSQ=(EOI/E0*(2.D0+EOI/E0))*ACON**2
         RG=sqrt(PSQ-PRG**2)/GAM ! R=(Ptheta)/gam
      ENDIF
      write(6,'('' Rguess,PRguess='',2f12.5)')RG,PRG
      FOE=0.0D+0
      bad=.false.

      if(ieo_type .eq. 2) then
         ieo_stp=nint(theo/(stp*tcon))
         if(allocated(reo)) deallocate(reo)
         allocate(reo(2,1,neo),stat=ierr)
         call allocate_check(ierr,'EO Storage reo')
      elseif(ieo_type .eq. 3) then
         ieo_stp=-1 ! inhibit part2 EO calcs
         allocate(cross(4*ndee),stat=ierr)
         call allocate_check(ierr,'EO Storage -cross ')
         if(allocated(EO_NUMBER).and.(size(EO_NUMBER).ne.4*ndee))&
              & deallocate(EO_NUMBER)
         if(.not. allocated(EO_NUMBER)) allocate (EO_NUMBER(4*ndee),stat&
              & =ierr)
         call allocate_check(ierr,'Gap Array EO_NUMBER')
! assume normal sprgapz use - 4 lines/dee then repeat ndee times
         eo_ngap=4 ! spiral gap behaviour
         do i=1,ndee
            k=4*(i-1)
            EO_NUMBER(k+1)=1
            EO_NUMBER(k+2)=2
            EO_NUMBER(k+3)=3
            EO_NUMBER(k+4)=4
         enddo
         if(allocated(reo)) deallocate(reo)
         allocate(reo(2,eo_ngap,neo),stat=ierr)
         call allocate_check(ierr,'EO Storage reo')
      elseif(ieo_type .eq. 4) then
         kstep_eo=nstp
         if(allocated(reop)) deallocate(reop)
         if(allocated(preop)) deallocate(preop)
         allocate(reop(kstep_eo,neo),stat=ierr)
         call allocate_check(ierr,'EO Storage reop')
         allocate(preop(kstep_eo,neo),stat=ierr)
         call allocate_check(ierr,'EO Storage reop')
         if(allocated(reo)) deallocate(reo)
         allocate(reo(2,1,4),stat=ierr)
         call allocate_check(ierr,'EO Storage reo')
      endif
      do IEO=1,NEO
         E=EOI+(IEO-1)*DEO
         PSQ=E/E0
         GAM=1.D+0+PSQ
         PSQ=PSQ*(2.D+0+PSQ)*ACON**2
         LastIter=.false.
         do ntry=1,20
         Y(1)=RG
         Y(2)=PRG
         Y(3)=1.D+0
         Y(4)=0.D+0
         Y(5)=0.D+0
         Y(6)=1.D+0
         Y(7)=0.D+0
         Y(8)=1.D+0
         Y(9)=0.D+0
         Y(10)=0.D+0
         Y(11)=1.D+0
         Y(12)=0.D+0
         NXCROS=0
         NZCROS=0
         Q=0.D+0
         do istp=1,nstp
            RST(ISTP)=y(1)
            PRST(ISTP)=y(2)
            yp5=y(5)
            yp10=y(10)
            do j=1,4
               ITH=2*(ISTP-1)+int(j/2)+1 !zero degs is ith=1
               IF (.not. BFLDCAL(ith,y(1),.true.,1) )then
                  bad=.true.
                  goto 10
               endif
               PTS=PSQ-Y(2)**2
               IF(PTS .lt. 0) Then
                  write(26,'(''PR=P'')')
                  bad=.true.
                  goto 10
               endif
               PT=SQRT(PTS)
               XC0=STP*Y(1)/PT
               CK(1)=XC0*Y(2)
               CK(2)=STP*(PT-Y(1)*BZ)
               XC1=STP*Y(2)/PT
               XC2=PSQ*XC0/PTS
               XC3=STP*(Y(1)*dbzdr+BZ)
               CK(3)=XC1*Y(3)+XC2*Y(4)
               CK(4)=-XC1*Y(4)-XC3*Y(3)
               CK(5)=XC1*Y(5)+XC2*Y(6)
               CK(6)=-XC1*Y(6)-XC3*Y(5)
 !              write(67,'(i5,f10.3,3e16.5)')ith,y(1),bz,dbzdr
               NEQ=6
               IF(lastIter) then ! PHASE, Z, AND RAV EQUATIONS (LAST TIME ONLY)
                  NEQ=12
                  XC4=STP*Y(1)*dbzdr-XC1*dbzdt
                  CK(7)=XC0*GAM-STP
                  CK(8)=XC0*Y(9)
                  CK(9)=XC4*Y(8)
                  CK(10)=XC0*Y(11)
                  CK(11)=XC4*Y(10)
                  CK(12)=STP*Y(1)
               endif
               DO  I=1,NEQ
                  RRK=ARK(J)*(CK(I)-BRK(J)*Q(I))
                  Y(I)=Y(I)+RRK
                  Q(I)=Q(I)+3.D+0*RRK+CRK(J)*CK(I)
               enddo !I
            enddo !J
            if(istp .gt. 1) then
               IF(Y(5)*YP5.LT.0.D+0)NXCROS=NXCROS+1
               IF(Y(10)*YP10.LT.0.D+0)NZCROS=NZCROS+1
            endif
          enddo ! istp
          RST(NSTP+1)=RST(1)
          PRST(NSTP+1)=PRST(1)
          if(LastIter) then
               exit !ntry
          else
             EP1=Y(1)-RG
             EP2=Y(2)-PRG
             DEN=(1.D+0-Y(3))*(1.D+0-Y(6))-Y(4)*Y(5)
             RG=RG+(EP1*(1.D+0-Y(6))+EP2*Y(5))/DEN
             PRG=PRG+(EP1*Y(4)+EP2*(1.D+0-Y(3)))/DEN
             IF(ABS(EP1)+ABS(EP2) .LT. 1.D-5*RG) LastIter=.true.
          endif
          enddo ! ntry
          if(ntry .ge. 20) then
             write(26,'(11H LOST   R =,F7.3,7H   PR =,F7.3)')rg,prg
             bad=.true.
             goto 10
          endif
          PHSP=PHS
          PHS=Y(7)/(SECT)
          RAV=Y(12)/(SECT)
          DO I=1,2
             I1=5*(I-1)
             CN=.5D+0*(Y(3+I1)+Y(6+I1))
             SN=1.-CN*CN
             if(SN .LT. 0) then
                RZNT=-ALOG(ABS(CN)+SQRT(-SN))
             else
                SN=SQRT(SN)
                IF(Y(5+I1) .LT. 0) SN=-SN
                RZNT=ATAN2(SN,CN)
                if(RZNT .lt. 0) RZNT=RZNT+6.283185
             endif
             RZNU(I)=RZNT/(SECT)
             arg(i)=rznt
          enddo ! i
          IF((E-EOI)/DEO.GE..5) FOE=FOE+3141.5927*DEO*(PHSP+PHS)
          if(ieo_type .gt. 1) then
            write(26,309,iostat=ierr)E,RAV,RG,PRG,PHS,RZNU(1),NXCROS,RZNU(2),NZCROS,FOE
309         FORMAT(4F10.5,F10.6,F10.5,I2,F10.5,I2,F10.3)
          endif
! now find gap crossings
          if(ieo_type .eq. 2) then ! fixed theta
             reo(1,1,ieo)=rst(ieo_stp+1)
             reo(2,1,ieo)=prst(ieo_stp+1)
          elseif(ieo_type .eq. 3) then ! spiral gaps
             do igp=1,eo_ngap
                found=.false.
                do i=0,nsc-1
                   do istp=1,nstp
                      k=istp-1+i*nstp
                      call thin(k,stp,rst(istp),prst(istp), &
                        rst(istp+1),prst(istp+1),psq,igp,fth,dummy,*10, &
                        cross(igp),reo(1,igp,ieo),reo(2,igp,ieo))
                      if(fth .ge. 0.0) then
                         found=.true.
                         exit
                      endif
                   enddo ! istp
                   if(found) then
                      exit
                   endif
                enddo ! i=1,nsc
                if(.not. found) then
                   write(6,'('' Missed Gap '',2i5)')igp,ieo
                   stop 'EO Error'
                endif
             enddo ! igp=1,eo_ngap
             write(26,313)(cross(J),(reo(i,j,ieo),i=1,2),J=1,eo_ngap)
313          FORMAT(1X,F10.3,2F10.5,F10.3,2F10.5,F10.3,2F10.5,1F10.3&
                  & ,2F10.5)
          elseif(ieo_type .eq. 4) then !all angles
             do i=1,kstep_eo
                REOP(i,ieo)=rst(i)
                PREOP(i,ieo)=prst(i)
             enddo
             write(26,313)(i*stp*tcon,rst(i),prst(i),i=1,kstep_eo)
          endif
      ENDDO !IEO

10    continue
      if(bad) then
         neo=ieo-1
         write(6,'("Last EO attempt failed")')
      endif
      write(6,'('' Number of EOs computed is '',2i5)')neo

!     EO DATA.  NOTE THAT R**2 AND PR/P ARE STORED
      if(ieo_type .eq. 1) then
        ! this is used to find eo or eigen ellipse
      else if(ieo_type .le. 2) then
         DO I=1,NEO
          PEO=(EOI+(I-1)*DEO)/E0
          PEO=ACON*SQRT(PEO*(2.+PEO))
          REO(1,1,I)=REO(1,1,I)**2
          REO(2,1,I)=REO(2,1,I)/PEO
        enddo
      else if(ieo_type .eq. 3) then
         DO I=1,NEO
          PEO=(EOI+(I-1)*DEO)/E0
          PEO=ACON*SQRT(PEO*(2.+PEO))
          do k=1,4
             REO(1,k,I)=REO(1,k,I)**2
             REO(2,k,I)=REO(2,k,I)/PEO
          enddo
        end do
      elseif(ieo_type .eq. 4) then
        do i=1,neo
            PEO=(EOI+(I-1)*DEO)/E0
            PEO=ACON*SQRT(PEO*(2.+PEO))
            do j=1,kstep_eo
              REOP(j,i)=REOP(j,I)**2
              PREOP(j,i)=PREOP(j,I)/PEO
            end do
        end do
      endif

      if(allocated(rst)) deallocate(rst,stat=ierr)
      if(allocated(prst)) deallocate(prst,stat=ierr)
      if(allocated(cross)) deallocate(cross,stat=ierr)
      return

      end subroutine eo_make

!************************************************
!
! called by parameter_load to generate the ellipse based starting conditions    
!
    subroutine fill_ellipse(r,pr,energy)
      real,intent(in) :: r ! rguess
      real,intent(in) :: pr ! rguess
      real,intent(in) :: energy ! energy to calculate EO     

      real :: area(2) ! ellipse area
      real :: a(2),b(2),co(2),si(2) ! ellipse shape, used when an eigen ellipse is generated
      real :: tt,bg ! for ellipse area calc
      integer :: neo_save,ity_save ! DUMMIES FOR STORING TEMP VALUES
      real :: eoi_save ! DUMMIES FOR STORING TEMP VALUES
      real*8:: gam ! relativistic gamma
      real:: rho,ratio,den ! for gaussian radius calc
      integer :: itype,i,ierr,j,n,l
      real, allocatable :: sigma(:,:),offsets(:)
      real, allocatable :: temp(:,:)
      integer :: iEllipseDim ! either 2 or 4
      integer :: nSave     
      
       TT=ACON*SQRT(Energy/E0*(2.+Energy/E0))*0.001/(conv*10)
       pScale=ACON*SQRT(Energy/E0*(2.+Energy/E0))*0.001
       GAM=Energy/E0+1.D0
       bg=sqrt(gam*gam-1.d+0) !beta*gamma
       if(UseEigen) then
          neo_save=neo
          eoi_save=eoi
          ity_save=ieo_type
          neo=1
          eoi=energy
          ieo_type=1
          if(deo .le. 0.0) deo=1.0
          call eo_make(r,pr)
          neo=neo_save
          eoi=eoi_save
          ieo_type=ity_save               
          if(ellipse_norm) then
            write(6,'('' The normalized emittances are '',2f10.5,&
                   ''pi mm-mrad'')')ellipse_area
            write(6,'('' The unnormalized emittances are '',2f10.5,&
                   ''pi mm-mrad'')')ellipse_area/bg
            area=ellipse_area/bg
          else
            write(6,'('' The normalized emittances are '',2f10.5,&
                   ''pi mm-mrad'')')ellipse_area*bg
            write(6,'('' The unnormalized emittances are '',2f10.5,&
                   ''pi mm-mrad'')')ellipse_area
            area=ellipse_area
          endif
          area=area*tt ! was entered as mm-mrad not length-length
          ! The ellipse centroid is the equilibrium orbit unless the UseEllipseOffset Flag is set
          if (.Not. UseEllipseOffset) then
            Ellipse_Offset_R = Y(1)
            Ellipse_Offset_Pr = Y(2)
          endif
          call get_AB(a,b,co,si,area)
       else
          write(6,'('' The normalized emittances are '',2f10.5,&
                 ''pi mm-mrad'')')ellipse_area*bg/tt
          write(6,'('' The unnormalized emittances are '',2f10.5,&
                 ''pi mm-mrad'')')ellipse_area/tt
          area=ellipse_area
          call Get_ABfromSigma(a,b,co,si,area,ellipse_type)
      endif
      ellipse_counter=1 !reset the counter used to run ellipse particles
      itype=ellipse_type
      select case (nEllipseFillType)
      case (1) ! edge
        if (itype < 3)then
          if(allocated(ellipse)) deallocate(ellipse)
          allocate(ellipse(2,nEllipseParticles),stat=ierr)
          call allocate_check(ierr,'Ellipse Storage')
          call FillEdge(a,b,si,co,Ellipse_Offset_R,Ellipse_Offset_Pr,itype)
          iEllipseDim=2
        else
          iEllipseDim=4
          n=nEllipseParticles-1.0
          if(n .lt. 8) n=8 ! fill edge looks for n<=9
          nSave=n*n+1
          n=n+1
          nEllipseParticles=n
          if(allocated(ellipse)) deallocate(ellipse)
          allocate(ellipse(4,nSave),stat=ierr)
          call allocate_check(ierr,'Ellipse Storage')
          allocate(temp(4,n),stat=ierr)
          call allocate_check(ierr,'Ellipse Temp Storage')
          call FillEdge(a,b,si,co,Ellipse_Offset_R,Ellipse_Offset_Pr,1)
          Temp(1:2,1:n-1)=ellipse(1:2,2:n)
          call FillEdge(a,b,si,co,Ellipse_Offset_R,Ellipse_Offset_Pr,2)
          Temp(3:4,1:n-1)=ellipse(1:2,2:n)
          ellipse(3:4,1)=0.0 ! first element is on axis
          l=1 ! skip the first element as in is the centre
          n=n-1
          do i=1,n
            do j=1,n
              l=l+1
              Ellipse(1:2,l)=Temp(1:2,i)
              Ellipse(3:4,l)=Temp(3:4,j)
            enddo
          enddo
          nEllipseParticles=nSave
          deallocate(temp)
        endif
      case (2) ! Gaussian
        if(GaussTruncate) then
          n=nEllipseParticles
          nEllipseParticles=nEllipseParticles/(1-2*erfc(sqrt(GaussianBound)))**2
          write(6,'(" truncating ellipse at ", f5.2," sigma, generating",i6,f8.5)')GaussianBound,nEllipseParticles,erfc(sqrt(GaussianBound))
        endif
        call ManageArrays(sigma,offsets,iEllipseDim,itype)
        call Get_Sigma(sigma,y,arg,area,itype,iEllipseDim)
        rho=sigma(1,2)/sqrt(sigma(1,1)*sigma(2,2))
        den=1/(1-rho**2)
        rho=2*sigma(1,2)/(sigma(1,1)*sigma(2,2))
        call FillGaussian(sigma,offsets,iEllipseDim)
        if(GaussTruncate) then
          l=0
          do i=1,nEllipseParticles
            ratio=(ellipse(1,i)-offsets(1))**2/sigma(1,1)+(ellipse(2,i)-offsets(2))**2/sigma(2,2)
            ratio=ratio-rho*(ellipse(1,i)-offsets(1))*(ellipse(2,i)-offsets(2))
            ratio=ratio*den
            !write(6,'(5f10.5)')ellipse(1,i),ellipse(2,i),ratio,den,rho
            if(ratio .le. GaussianBound) then
              l=l+1
              ellipse(1:2,l)=ellipse(1:2,i)
              if(l .ge. n) then
                exit
              endif
            endif
          enddo
          nEllipseParticles=l
        endif
        write(6,'(" Number of remaining particles is ",i6)')nEllipseParticles
      case (3) ! uniform
        if(itype .eq. 3) then
          n=nEllipseParticles
          nSave=n*n+1
          !n=n+1
          nEllipseParticles=nsave ! set to create the correct array size
          call ManageArrays(sigma,offsets,iEllipseDim,itype)
          call Get_Sigma(sigma,y,arg,area,itype,iEllipseDim)          
          write(6,'(" Uniform Sigma")')
          write(6,'(4f10.4)')sigma
          write(6,'(" Offsets ",4f10.4)')offsets
          write(6,*)
          nEllipseParticles=n !set to generate the correct number in each plane
          call FillUniform(sigma(1:2,1:2),offsets(1:2),ellipse,0)
          allocate(temp(4,n),stat=ierr)
          call allocate_check(ierr,'Ellipse Temp Storage')
          Temp(1:2,1:n)=ellipse(1:2,1:n)
          call FillUniform(sigma(3:4,3:4),offsets(3:4),ellipse,0)
          Temp(3:4,1:n)=ellipse(1:2,1:n)
          ellipse(3:4,1)=0.0 ! first element is on axis
          l=0 ! skip the first element as in is the centre
          n=nEllipseParticles ! this was reset by FillUniform
          do i=1,n
            do j=1,n
              l=l+1
              Ellipse(1:2,l)=Temp(1:2,i)
              Ellipse(3:4,l)=Temp(3:4,j)
            enddo
          enddo
          nEllipseParticles=n*n
          deallocate(temp)
        else
          call ManageArrays(sigma,offsets,iEllipseDim,itype)
          write(6,'(" Uniform Sigma")')
          write(6,'(2f10.4)')sigma
          write(6,'(" Offsets ",2f10.4)')offsets
          write(6,*)
          call FillUniform(sigma,offsets,ellipse,0)
        endif
      case default
          stop "Illegal Ellipse Fill"
      end select
      ierr=openFile(52,return_msg)
      if (varyTau) then
        call FillTau()
        do i=1,nEllipseParticles
          write(52,'(5f12.5)')(ellipse(j,i),ellipse(J+1,i),j=1,iEllipseDim,2),tau(i)
        end do
      else
        do i=1,nEllipseParticles
          write(52,'(4f12.5)')(ellipse(j,i),ellipse(J+1,i),j=1,iEllipseDim,2)
        end do        
      endif
      if(allocated(sigma)) deallocate(sigma)
      if(allocated(offsets)) deallocate(offsets)
    end subroutine
!************************************************
!
    subroutine ManageArrays(sigma,offsets,iEllipseDim,itype)
      real, allocatable,intent(inout) :: sigma(:,:),offsets(:)
      integer,intent(out):: iEllipseDim
      integer,intent(in):: itype ! ellipse type x,z, or both
      ! local values
      integer:: ierr

      if(allocated(sigma)) deallocate(sigma)
      if(allocated(offsets)) deallocate(offsets)
      if(allocated(ellipse)) deallocate(ellipse)
              
      select case (itype)
      case (3)
        allocate(ellipse(4,nEllipseParticles),stat=ierr)
        call allocate_check(ierr,'Ellipse Storage')
        ellipse=0.0            
        allocate(sigma(4,4),stat=ierr)
        call allocate_check(ierr,'Ellipse Sigma Storage')
        allocate(offsets(4),stat=ierr)
        call allocate_check(ierr,'Ellipse Offset Storage')
        sigma=0.0
        offsets(1)=Ellipse_Offset_R
        offsets(2)=Ellipse_Offset_Pr
        offsets(3)=Ellipse_Offset_Z
        offsets(4)=Ellipse_Offset_Pz
        iEllipseDim=4
      case (1,2)
        allocate(ellipse(2,nEllipseParticles),stat=ierr)
        call allocate_check(ierr,'Ellipse Storage')
        allocate(sigma(2,2),stat=ierr)
        call allocate_check(ierr,'Ellipse Sigma Storage')
        allocate(offsets(2),stat=ierr)
        call allocate_check(ierr,'Ellipse Offset Storage')
        ellipse=0.0
        sigma=0
        offsets=0
        if (itype .eq. 1) then
            offsets(1)=Ellipse_Offset_R
            offsets(2)=Ellipse_Offset_Pr
        else
            offsets(1)=Ellipse_Offset_Z
            offsets(2)=Ellipse_Offset_Pz
        endif
        iEllipseDim=2
      end select
    end subroutine ManageArrays

!************************************************
!
    SUBROUTINE Get_ABfromSigma(xa,xb,xco,xsi,area,itype)

! Dummy Variables
        real, intent(out) ::xa(2),xb(2),xco(2),xsi(2)
        real,intent (in) :: area(2) ! area is assume to be in length-length units as this is set before call to eigen
        integer,intent(in) :: itype

! sigma values assumed to be in mm-mrad
        
! local variables
        real ::al,be,ga,con,de,s
        real :: eps
        real ::a,b,co,si

        if(itype .eq. 3) then
          !write(6,'("area,s11,s22,s12",5f10.5)')area(1),sigma_in(1),sigma_in(2),sigma_in(3)
          ga = sigma_in(2)/area(1)
          be = sigma_in(1)/area(1)
          if (sigma_in(3) .ne. 0) then
            al = -sigma_in(3)/area(1)
          else
            al=0.0
          endif
          CON=SQRT((BE-GA)**2+4.*AL*AL)
          B=(SQRT(CON*CON+4.)-CON)/2.
          B=SQRT(B*AREA(1))
          DE=1.
          IF(AL.LT.0.) DE=-1.
          if (CON .ne. 0.0) Then
            SI=DE*SQRT((1.+(BE-GA)/CON)/2.)
            CO=SQRT(1.-SI*SI)
          else
            SI=0.0
            CO=1.0
          endif
          A=AREA(1)/B
          !write(6,'("alpha,beta,gamma,a,b",6f12.5)')al,be,ga,a,b
          xa(1)=a
          xb(1)=b
          xco(1)=co
          xsi(1)=si
          
          !write(6,'("area,s11,s22,s12",5f10.5)')area(2),sigma_in(4),sigma_in(5),sigma_in(6)
          ga = sigma_in(5)/area(2)
          be = sigma_in(4)/area(2)
          if (sigma_in(6) .ne. 0) then
            al = -sigma_in(6)/area(2)
          else
            al=0.0
          endif
          CON=SQRT((BE-GA)**2+4.*AL*AL)
          B=(SQRT(CON*CON+4.)-CON)/2.
          B=SQRT(B*AREA(2))
          DE=1.
          IF(AL.LT.0.) DE=-1.
          if (CON .ne. 0.0) Then
            SI=DE*SQRT((1.+(BE-GA)/CON)/2.)
            CO=SQRT(1.-SI*SI)
          else
            SI=0.0
            CO=1.0
          endif
          A=AREA(2)/B
          !write(6,'("alpha,beta,gamma,a,b",6f12.5)')al,be,ga,a,b
          xa(2)=a
          xb(2)=b
          xco(2)=co
          xsi(2)=si          
        else
          !write(6,'("area,s11,s22,s12",5f10.5)')area(itype),sigma_in(1),sigma_in(2),sigma_in(3)
          ga = sigma_in(2)/area(itype)
          be = sigma_in(1)/area(itype)
          if (sigma_in(3) .ne. 0) then
            al = -sigma_in(3)/area(itype)
          else
            al=0.0
          endif
          CON=SQRT((BE-GA)**2+4.*AL*AL)
          B=(SQRT(CON*CON+4.)-CON)/2.
          B=SQRT(B*AREA(itype))
          DE=1.
          IF(AL.LT.0.) DE=-1.
          if (CON .ne. 0.0) Then
            SI=DE*SQRT((1.+(BE-GA)/CON)/2.)
            CO=SQRT(1.-SI*SI)
          else
            SI=0.0
            CO=1.0
          endif
          A=AREA(itype)/B
          !write(6,'("alpha,beta,gamma,a,b",6f12.5)')al,be,ga,a,b
          xa(1)=a
          xb(1)=b
          xco(1)=co
          xsi(1)=si
        endif

    end subroutine
!************************************************
!
    SUBROUTINE Get_AB(a,b,co,si,area)

! Dummy Variables
        real,intent(in) :: area(2)
        real, intent(out) ::a(2),b(2),si(2),co(2)

! local variables
        real ::al,be,ga,con,de,s

        ! itype=1 is x
          call get_abs(a(1),b(1),co(1),si(1),area,1)
        ! itype=2 is z
          call get_abs(a(2),b(2),co(2),si(2),area,2)
        return
    end SUBROUTINE Get_AB
        
!************************************************
!
    SUBROUTINE Get_ABS(a,b,co,si,area,itype)

! Dummy Variables
        real,intent(in) :: area(2)
        real, intent(out) ::a,b,si,co
        integer,intent(in) :: itype

! local variables
        real ::al,be,ga,con,de,s

        S=SIN(ARG(Itype))
        AL=(Y(3+5*(Itype-1))-Y(6+5*(Itype-1)))/2./S
        BE=Y(5+5*(Itype-1))/S
        GA=(1.+AL*AL)/BE
        CON=SQRT((BE-GA)**2+4.*AL*AL)
        B=(SQRT(CON*CON+4.)-CON)/2.
        B=SQRT(B*AREA(itype))
        DE=1.
        IF(AL.LT.0.) DE=-1.
        SI=DE*SQRT((1.+(BE-GA)/CON)/2.)
        CO=SQRT(1.-SI*SI)
        A=AREA(itype)/B
        return
    end SUBROUTINE Get_ABS
!
!********************************************************
!
    SUBROUTINE Get_Sigma(sigma,y,arg,area,itype,ndim)

! Dummy Variables
        integer,intent(in) ::ndim
        real*8,intent(in) :: Y(12)
        real,intent(in) :: arg(2)
        real,intent(in) :: area(2)
        integer,intent(in) :: itype
        real, intent(out) :: sigma(ndim,ndim)

! local variables
        real ::al,be,ga,con,de,s
        real ::a,b,si,co
        
        sigma=0.0
        if (UseEigen) then
            if (itype .eq. 3)then ! both
                S=SIN(ARG(1))
                AL=(Y(3)-Y(6))/2./S
                BE=Y(5)/S
                GA=(1.+AL*AL)/BE
                sigma(1,1) = area(1)/ga
                sigma(2,2) = area(1)/be
                sigma(1,2) = area(1)*al
                sigma(2,1) = sigma(1,2)
                ! Now z
                S=SIN(ARG(2))
                AL=(Y(3+5)-Y(6+5))/2./S
                BE=Y(5+5)/S
                GA=(1.+AL*AL)/BE
                sigma(3,3) = area(2)/ga
                sigma(4,4) = area(2)/be
                sigma(3,4) = area(2)*al
                sigma(4,3) = sigma(3,4)
            else
                S=SIN(ARG(Itype))
                AL=(Y(3+5*(Itype-1))-Y(6+5*(Itype-1)))/2./S
                BE=Y(5+5*(Itype-1))/S
                GA=(1.+AL*AL)/BE
                sigma(1,1) = area(itype)/ga
                sigma(2,2) = area(itype)/be
                sigma(1,2) = area(itype)*al
                sigma(1,2) = sigma(2,1)
            endif

        else ! This means the sigma values were entered by the user
            if(itype .eq. 3) then ! both x & z
                sigma(1,1)=sigma_in(1)
                sigma(2,2)=sigma_in(2)
                sigma(1,2)=sigma_in(3)
                sigma(2,1)=sigma_in(3)
                sigma(3,3)=sigma_in(4)
                sigma(4,4)=sigma_in(5)
                sigma(4,3)=sigma_in(6)
                sigma(3,4)=sigma_in(6)
                if(UseCorr) then
                  sigma(1,3)=sigma_in(7)
                  sigma(3,1)=sigma_in(7)
                  sigma(1,4)=sigma_in(8)
                  sigma(4,1)=sigma_in(8)
                  sigma(2,3)=sigma_in(9)
                  sigma(3,2)=sigma_in(9)
                  sigma(2,4)=sigma_in(10)
                  sigma(4,2)=sigma_in(10)
                endif
            else ! either x or z
                sigma(1,1)=sigma_in(1)
                sigma(2,2)=sigma_in(2)
                sigma(1,2)=sigma_in(3)
                sigma(2,1)=sigma_in(3)
            endif
        endif
        return
    end SUBROUTINE Get_Sigma

!************************************************
!
! called by parmameter_load to get EO based starting conditions
!
    subroutine get_EO(r,pr,energy)
      real,intent(in) :: r ! rguess
      real,intent(in) :: pr ! rguess
      real,intent(in) :: energy ! energy to calculate EO
      
      integer :: neo_save,ity_save
      real :: eoi_save
      
      integer:: ierr ! error code from array allocation
            
        neo_save=neo
        eoi_save=eoi
        ity_save=ieo_type
        neo=1
        eoi=energy
        ieo_type=1
        if(deo .le. 0.0) deo=1.0
       
        call eo_make(r,pr)
        
        if(allocated(ellipse)) deallocate(ellipse)
        allocate(ellipse(2,1),stat=ierr)
        call allocate_check(ierr,'EO Storage reo')
        ! at this point the offsets are x, px
        if (UseEllipseOffset) then
          ellipse(1,1)=Y(1) + ellipse_offset_r
          ellipse(2,1)=Y(2) + ellipse_offset_pr
        else
          ellipse(1,1)=Y(1)
          ellipse(2,1)=Y(2)
        endif
        ! to maintain the normal condition these are now r,pr for ellipse if required
        call setEllipseOffset(ellipse(1,1),ellipse(2,1),ellipse_offset_z,ellipse_offset_pz,.true.)
        !reset the counter used to run ellipse particles
        ellipse_counter=1 
        ! return orginal settings for EO calcs
        neo=neo_save
        eoi=eoi_save
        ieo_type=ity_save
    end subroutine
!************************************************
      subroutine load_ellipse_info(rin,prin,zin,pzin,tin,forceTheta)
! dummy variables
      real  ,intent(inout) :: rin,prin,zin,pzin,tin
      logical, intent(out) :: forceTheta

      if(ellipse_type .eq. 0) then
         rin=ellipse(1,1)
         prin=ellipse(2,1)
         run_ellipse=.false.
      elseif(ellipse_type .eq. 1) then
         rin=ellipse(1,ellipse_counter)
         prin=ellipse(2,ellipse_counter)
      elseif(ellipse_type .eq. 2) then
         if(ellipse_counter .eq. 1) ellipse_counter=2
         zin=ellipse(1,ellipse_counter)
         pzin=ellipse(2,ellipse_counter)
         rin=Ellipse_Offset_R
         prin=Ellipse_Offset_Pr
      elseif(ellipse_type .eq. 3) then
         rin=ellipse(1,ellipse_counter)
         prin=ellipse(2,ellipse_counter)
         zin=ellipse(3,ellipse_counter)
         pzin=ellipse(4,ellipse_counter)
      endif
      if (VaryTau) then
        tin = tau(ellipse_counter)
      endif
      if(debug1) write(6,'("Loaded ellipse particle # ",i5," r= ",f12.5)')ellipse_counter,rin
      ellipse_counter=ellipse_counter+1
      if(ellipse_counter .eq. nEllipseParticles+1) run_ellipse=.false.
      forceTheta = UseEigen ! at present eigen ellipse only available at th=0

      end subroutine load_ellipse_info
      
!************************************************
!  Generate tau values for the ellipse particles - uniform distribution
!
    Subroutine FillTau()
! Dummy Variables

! Local Variables
        integer :: ierr
        
        if(allocated(tau)) deallocate(tau)
        allocate(tau(nEllipseParticles),stat=ierr)
        call allocate_check(ierr,'Tau Storage')
        CALL RANDOM_NUMBER (tau)
        tau = TauMin+TauRange*tau
    
    end subroutine FillTau
    
!************************************************
!  Generate particles on the edge of ellipse
!
    Subroutine FillEdge(a,b,si,co,r,pr,itype)
! Dummy Variables
        real, intent(in) :: a(2),b(2),si(2),co(2)
        real*8, intent(in) :: r,pr
        integer, intent(in) :: itype
! local variables      
        real ::x(8),p(8), dth,th
        integer :: i
 
     write(6,'(" a,b,si,co ",4f10.5)')a(itype),b(itype),si(itype),co(itype)
        X(1)=-A(itype)*SI(itype)
        P(1)=A(itype)*CO(itype)
        X(5) =-X(1)
        P(5)=-P(1)
        X(3)=B(itype)*CO(itype)
        P(3)=B(itype)*SI(itype)
        X(7)=-X(3)
        P(7)=-P(3)
        X(2)=(X(1)+X(3))/1.41421
        P(2)=(P(1)+P(3))/1.41421
        X(6)=-X(2)
        P(6)=-P(2)
        X(8)=(X(1)-X(3))/1.41421
        P(8)=(P(1)-P(3))/1.41421
        X(4)=-X(8)
        P(4)=-P(8)

        !write(6,'('' R0,PR0=''2f12.5)')r,pr
        if(itype .eq. 2) then
          ellipse(1,1)=0.0
          ellipse(2,1)=0.0       
        else
          ellipse(1,1)=r
          ellipse(2,1)=pr
        endif
        if (nEllipseParticles .le. 9) then
            do i=1,8
               if(itype .eq. 1) then
                  ellipse(1,i+1)=x(i)+r
                  ellipse(2,i+1)=p(i)+pr
                  !write(6,'(" R=",f12.2," PR =",f12.5)')ellipse(1,i+1),ellipse(2,i+1)
               else
                  ellipse(1,i+1)=x(i)
                  ellipse(2,i+1)=p(i)
                  !write(6,'(" Z=",f12.5," PZ =",f12.5)')ellipse(1,i+1),ellipse(2,i+1)
               endif
            enddo
        else
            dth = tpi/(nEllipseParticles-1)
            do i=1,nEllipseParticles-1
                th = dth*(i-1) ! angle around ellipse
                if(itype .eq. 1) then !x
                  ellipse(1,i+1)=r  -a(itype)*si(itype)*cos(th) + b(itype)*co(itype)*sin(th)
                  ellipse(2,i+1)=pr +a(itype)*co(itype)*cos(th) + b(itype)*si(itype)*sin(th)
                  !write(6,'('' R,PR ='',2f12.5)')ellipse(1,i+1),ellipse(2,i+1)
                else !z
                  ellipse(1,i+1)= -a(itype)*si(itype)*cos(th) + b(itype)*co(itype)*sin(th)
                  ellipse(2,i+1)= +a(itype)*co(itype)*cos(th) + b(itype)*si(itype)*sin(th)
                  !write(6,'('' Z,PZ ='',2f12.5)')ellipse(1,i+1),ellipse(2,i+1)
                endif
            end do
         endif
        return
    end subroutine FillEdge

!************************************************
!  Generate uniform distribution of particles in ellipse
!
    subroutine FillUniform(sigma,offsets,ellipse,ioff)
      real,intent(in) :: sigma(2,2) ! The sigma matrix for the ellipse
      real,intent(in) :: offsets(2) ! Centre of the ellipse
      real,intent(out) :: ellipse(4,*)
      integer,intent(in) :: ioff ! offset in ellipse
            
      real:: r1 ! radial step size in unit circle
      integer:: i ! number of radial steps
      real,allocatable:: rr(:),th(:)
      integer:: k ! counter
      integer::j,l,ll,ierr
      real:: m11,m12,m21,m22
      real:: eps
      real:: xx,yy
      
      r1=sqrt(1.0/nEllipseParticles)
      i=sqrt(float(nEllipseParticles))

      allocate(rr(nEllipseParticles),stat=ierr)
      call allocate_check(ierr,'Ellipse r unit storage')
      allocate(th(nEllipseParticles),stat=ierr)
      call allocate_check(ierr,'Ellipse th unit storage')

! generate equal spacing on a unit circle      
      k=1
      rr(1)=0.
      th(1)=0.
      do l=1,i-1
        ll=2*l+1
        write(99,'(3i5)')l,ll,k
        do j=1,ll
          k=k+1
          rr(k)=(l+0.5)*r1
          th(k)=(j-1)*tpi/ll
        enddo !j
      enddo !l
      nEllipseParticles=k

      eps=sqrt(Sigma(1,1)*Sigma(2,2)-Sigma(1,2)**2)
      M11=sqrt(Sigma(1,1)/eps)
      M12=0
      M21=Sigma(1,2)/sqrt(Sigma(1,1)*eps)
      M22=1/M11
      rr=sqrt(eps)*rr
      
      do j=1,nEllipseParticles
        ellipse(1+ioff,j)=M11*rr(j)*cos(th(j))+offsets(1)
        ellipse(2+ioff,j)=M21*rr(j)*cos(th(j))+M22*rr(j)*sin(th(j))+offsets(2)
      enddo
      if(allocated(rr)) deallocate(rr)
      if(allocated(th)) deallocate(th)
      return

    end subroutine FillUniform
    
!************************************************
!  Generate Gaussian distribution of particles in ellipse
!
    subroutine FillGaussian(c,a,ndim) !sigma,offset,dim

! Dummy variables
    integer, intent(in) ::ndim ! number of independant variables, normally 2 or 4
    Real(kind=4), intent(in) :: c(ndim,ndim) ! The sigma matrix for the ellipses
    Real(kind=4), intent(in) :: a(ndim) ! The centroid positions

! Local variables
      integer(kind=4) i,j,k
      integer(kind=4) errcode
      integer info
      integer brng,method
      Integer seed/7777777/
      integer me
      real(kind=4) t(ndim,ndim)
      real(kind=4) tt(ndim,ndim)

      write(6,'(" Gaussian Sigma")')
      write(6,'(<ndim>f10.4)')c
      write(6,'("Offsets ",<ndim>f10.4)')a
      write(6,*)
#if SIGMA == 1    
      brng=VSL_BRNG_MCG31
      !seed=7777777
      method=VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
      me=VSL_MATRIX_STORAGE_FULL
#endif
      t  = c

      if (debug6) then
!     Variance-covariance matrix for test (should be symmetric,positive-definite)

          print *,'Variance-covariance matrix C'
          write (*,'(<ndim>F7.3)') c
          print *,''

          print *,'Mean vector a:'
          write (*,'(<ndim>F7.3)') a
          print *,''
      endif
#if SIGMA == 1
      call spotrf('U',ndim,t,ndim,info)
#endif
     

!     **************************************************************************
! Clear the upper half of t
      do i=1,ndim
        do j=1,ndim
          if(i>j) t(i,j)=0.0
        end do
      end do

      if(debug6) then
          do i=1,ndim
            do j=1,ndim
              tt(j,i)=0.0
              do k=1,ndim
                tt(j,i)=tt(j,i)+t(k,i)*t(k,j)
              end do
            end do
          end do

          print *,'SPOTRF subroutine'
          print *,'*****************'
          print *,'Matrix T:             Matrix T*T (should be equal to C):'
          do i=1,ndim
            print 10,t(:,i),tt(:,i)
          end do
    10    format(<ndim>F7.3,'  ',<ndim>F7.3)
          print *,''
      endif
!     **************************************************************************

#if SIGMA == 1
!     Stream initialization
      if (.not. StreamInitialized) then
          errcode=vslnewstream(stream,brng,seed)
          call CheckVslError(errcode)
          StreamInitialized = .true.
      endif

!     Generating random numbers from multivariate normal distribution
      errcode=vsrnggaussianmv(method,stream,nEllipseParticles,ellipse,ndim,me,a,t)
      call CheckVslError(errcode)
#else
      write(6,'("************* Guassian routines not compiled *************")')
#endif

      return
    end subroutine FillGaussian

!************************************************
! Configure the ellipse information - called by parameter_load
    logical function SetEllipseType()
    
        logical :: error
        integer :: ieo_com
        character*10 ellipse_commands(7)/'XE','ZE','BOTHE','XI','ZI','BOTHI','CORR'/

        real,external :: rdata_in
        integer,external :: idata_in,cdata_in
        logical,external :: ldata_in
      
        error = .false.
        ieo_com=cdata_in(.false.,ellipse_commands,7,0)
        if(ieo_com .le. 0) then
          SetEllipseType =.true. ! in error
          return
        endif

        select case (ieo_com)
        case (1) ! x eigen
            ellipse_type = 1
            UseEigen = .true.
            ellipse_area(1)=rdata_in(.false.,1.e-10,1e6,0.0)
            ellipse_area(2) =0.0
            if(ellipse_area(1) .le. 0.) then
                write(6,'(" Ellipse area is zero ")')
                error=.true.
            endif
            ellipse_norm=ldata_in(.false.,.true.)
        case (2) ! z eigen
            ellipse_type =2
            UseEigen = .true.
            ellipse_area(2)=rdata_in(.false.,1.e-10,1e6,0.0)
            ellipse_area(1) =0.0
            if(ellipse_area(2) .le. 0.) then
                write(6,'(" Ellipse area is zero ")')
                error=.true.
            endif
            ellipse_norm=ldata_in(.false.,.true.)
        case (3) ! both eigen
            ellipse_type = 3
            UseEigen = .true.
            ellipse_area(1)=rdata_in(.false.,1.e-10,1e6,0.0)
            ellipse_area(2)=rdata_in(.false.,1.e-10,1e6,0.0)
            if(ellipse_area(1)*ellipse_area(2) .le. 0.) then
                write(6,'(" Ellipse area is zero ")')
                error=.true.
            endif
            ellipse_norm=ldata_in(.false.,.true.)
        case (4) ! x input
            ellipse_type = 1
            UseEigen = .false.
            sigma_in(1)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(2)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(3)=rdata_in(.false.,-1.e6,1e6,0.0)
            if(sigma_in(3)**2 .ge. sigma_in(1)*sigma_in(2)) then
              write(6,'(" Ellipse area zero ",3f12.5)')sigma_in(1),sigma_in(2),sigma_in(3)
              error=.true.
            endif
            ellipse_area(1)=sqrt(sigma_in(1)*sigma_in(2)-sigma_in(3)**2)
            ellipse_area(2)=0.0
            ellipse_norm=.false.
        case (5) ! z input
            ellipse_type =2
            UseEigen = .false.
            sigma_in(1)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(2)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(3)=rdata_in(.false.,-1.e6,1e6,0.0)
            if(sigma_in(3)**2 .ge. (sigma_in(1)*sigma_in(2)) ) then
              write(6,'(" Ellipse area zero ",3f12.5)')sigma_in(1),sigma_in(2),sigma_in(3)
              error=.true.
            endif
            ellipse_area(2)=sqrt(sigma_in(1)*sigma_in(2)-sigma_in(3)**2)
            ellipse_area(1)=0.0
            ellipse_norm=.false.
        case (6) ! both input
            ellipse_type = 3
            UseEigen = .false.
            sigma_in(1)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(2)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(3)=rdata_in(.false.,-1.e6,1e6,0.0)
            sigma_in(4)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(5)=rdata_in(.false.,1.e-10,1e6,0.0)
            sigma_in(6)=rdata_in(.false.,-1.e6,1e6,0.0)
            if(sigma_in(3)**2 .ge. sigma_in(1)*sigma_in(2)) then
              write(6,'(" x Ellipse area zero ",3f12.5)')sigma_in(1),sigma_in(2),sigma_in(3)
              error=.true.
            endif
            if(sigma_in(6)**2 .ge. sigma_in(4)*sigma_in(5)) then
              write(6,'(" z Ellipse area zero ",3f12.5)')sigma_in(4),sigma_in(5),sigma_in(6)
              error=.true.
            endif
            ellipse_area(1)=sqrt(sigma_in(1)*sigma_in(2)-sigma_in(3)**2)
            ellipse_area(2)=sqrt(sigma_in(4)*sigma_in(5)-sigma_in(6)**2)
            ellipse_norm=.false.
        case (7) ! corr
            if (ellipse_type .ne. 3) then
              write(6,'("Cannot add a correlation without BOTHI set")')
              !stop 'Bad input ellipse correlation'
              return_msg='Bad input ellipse correlation'
              error=.true.
              ierror=100
            endif
            sigma_in(7)=rdata_in(.false.,-1.e6,1e6,0.0)
            sigma_in(8)=rdata_in(.false.,-1.e6,1e6,0.0)
            sigma_in(9)=rdata_in(.false.,-1.e6,1e6,0.0)
            sigma_in(10)=rdata_in(.false.,-1.e6,1e6,0.0)
            UseCorr = .true.
        case default
           write(6,'('' Illegal Ellipse type'')')
           error=.true.
        end select
        SetEllipseType = error
    end function SetEllipseType
!************************************************
!  Set a truncation boundry for gaussian fill type
    logical function SetEllipseBound()
      real,external :: rdata_in
      
      GaussianBound=rdata_in(.false.,0.,100.,0.)
      if (GaussianBound .ne. 0) then
        GaussTruncate=.true.
      else
        GaussTruncate=.false.
      endif
    
      SetEllipseBound=.false.
    end function
!************************************************
!  Set the Fill Type
    logical function SetEllipseFillType()

      integer,external :: cdata_in
      character*10 ellipse_fill_com(3)/'EDGE','GAUSSIAN','UNIFORM'/
      
        nEllipseFillType = cdata_in(.false.,ellipse_fill_com,3,1)
        SetEllipseFillType = .false.
    end function SetEllipseFillType
    
!************************************************    
!  Set the number of particles in the ellipse
    logical Function SetEllipseNumber()

      integer,external :: idata_in

        nEllipseParticles = idata_in(.false.,1,100000,9)
        SetEllipseNumber = .false.
    End Function SetEllipseNumber
    
!************************************************    
!  Read the ellipse offset from the input stream
!  note: For type 0 the offset is x,px otherwise it is r,pr
    Logical Function ReadEllipseOffset
    
        real,external :: rdata_in
        
        Ellipse_Offset_R = rdata_in(.false.,-1.e6,1.e6,0.0)
        Ellipse_Offset_PR = rdata_in(.false.,-1.e6,1.e6,0.0)
        Ellipse_Offset_Z = rdata_in(.false.,-1.e6,1.e6,0.0)
        Ellipse_Offset_PZ = rdata_in(.false.,-1.e6,1.e6,0.0)
        UseEllipseOffSet = .true.
        write(6,'(" Read Ellipse offset ",4f10.5,l2)')Ellipse_Offset_R,Ellipse_Offset_PR,Ellipse_Offset_Z,Ellipse_Offset_PZ
        ReadEllipseOffset=.false.
    End Function ReadEllipseOffset  

!************************************************    
!  Set the ellipse offset
!  note: For type 0 the offset is x,px otherwise it is r,pr
    Subroutine SetEllipseOffset(r,pr,z,pz,offset)
    real, Intent(in) :: r,pr,z,pz
    logical , intent(in) ::offset
        Ellipse_Offset_R = r
        Ellipse_Offset_PR = pr
        Ellipse_Offset_Z = z
        Ellipse_Offset_PZ = pz
        UseEllipseOffSet = offset
        write(6,'(" Set Ellipse offset ",4f10.5,l2)')Ellipse_Offset_R,Ellipse_Offset_PR,Ellipse_Offset_Z,Ellipse_Offset_PZ
    End Subroutine SetEllipseOffset   

!************************************************
! Random Tau for the ellipse - called by parameter_load

    logical function SetVaryTau

        integer :: iv_com
        character*10 vary_commands(3)/'N','P','T'/

        real,external :: rdata_in
        integer,external :: idata_in,cdata_in
!       logical,external :: ldata_in
      
        iv_com=cdata_in(.false.,vary_commands,3,0)
        
        select case (iv_com)
          case (1) !none
            VaryTau = .false.
          case (2) !Phi
            VaryTau = .true.
            VaryPot=0.0
          case (3) !Tau
            VaryTau = .true.
            VaryPot=1.0
          case default
            write(6,'(" Illegal Vary Command")')
            SetVaryTau = .true.
        end select
        TauMin = rdata_in(.false.,-360.,360.,0.)
        TauRange = rdata_in(.false.,0.,360.,10.)
        SetVaryTau = .false.
    end function SetVaryTau

!************************************************
! Get the current setting or POT switch - called by parameter_load
! if not Varying then use default otherwise value is stored locally

  real function GetPot(Def)
    ! dummy variables
    real, intent(in) :: Def
    
    if (VaryTau) then
      GetPot=VaryPot
    else
      GetPOt=Def
    endif
    
  end function GetPot
  
  subroutine clean_eo_info

    if(allocated(reop)) deallocate(reop)
    if(allocated(preop)) deallocate(preop)
    if(allocated(REO)) deallocate(REO)
    if(allocated(EO_NUMBER)) deallocate(EO_NUMBER)
    if(allocated(ellipse)) deallocate(ellipse)
    if(allocated(tau)) deallocate(tau)
  end subroutine
    
!************************************************    
! Intel VSL Error codes from MKL library   
     SUBROUTINE CheckVslError( num )

        USE MKL_VSL_TYPE
        INTEGER(KIND=4) :: num

      if ( num == VSL_ERROR_CPU_NOT_SUPPORTED ) then
        print 33, "Error: CPU version is not supported (code ",          &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_FEATURE_NOT_IMPLEMENTED ) then
        print 33, "Error: this feature not implemented yet (code ",      &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_UNKNOWN ) then
        print 33, "Error: unknown error (code ",num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_BADARGS ) then
        print 33, "Error: bad arguments (code ",num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_MEM_FAILURE ) then
        print 33, "Error: memory failure. Memory allocation problem",    &
     &            " maybe (code ",num,")."
        stop 1
      end if

      if ( num == VSL_ERROR_NULL_PTR ) then
        print 33, "Error: null pointer (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_INVALID_BRNG_INDEX ) then
        print 33, "Error: invalid BRNG index (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED ) then
        print 33, "Error: leapfrog initialization is unsupported",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED ) then
        print 33, "Error: skipahead initialization is unsupported",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BRNGS_INCOMPATIBLE ) then
        print 33, "Error: BRNGs are not compatible for the operation",   &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_STREAM ) then
        print 33, "Error: random stream is invalid (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BRNG_TABLE_FULL ) then
        print 33, "Error: table of registered BRNGs is full (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE ) then
        print 33, "Error: value in StreamStateSize field is bad",        &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_WORD_SIZE ) then
        print 33, "Error: value in WordSize field is bad (code ",        &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_NSEEDS ) then
        print 33, "Error: value in NSeeds field is bad (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_NBITS ) then
        print 33, "Error: value in NBits field is bad (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_UPDATE ) then
        print 33, "Error: number of updated entries in buffer is",       &
     &            " invalid (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_NO_NUMBERS ) then
        print 33, "Error: zero number of updated entries in buffer",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM ) then
        print 33, "Error: abstract random stream is invalid (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_CLOSE ) then
        print 33, "Error: can`t close file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_OPEN ) then
        print 33, "Error: can`t open file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_WRITE ) then
        print 33, "Error: can`t write to file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_FILE_READ ) then
        print 33, "Error: can`t read from file (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_BAD_FILE_FORMAT ) then
        print 33, "Error: file format is unknown (code ",num,")."
        stop 1
      end if

      if ( num == VSL_RNG_ERROR_UNSUPPORTED_FILE_VER ) then
        print 33, "Error: unsupported file version (code ",num,")."
        stop 1
      end if

            if ( num == VSL_SS_ERROR_ALLOCATION_FAILURE ) then
        print 33, "Error: memory allocation failure (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_DIMEN ) then
        print 33, "Error: bad dimension value (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OBSERV_N ) then
        print 33, "Error: bad number of observations (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_STORAGE_NOT_SUPPORTED ) then
        print 33, "Error: storage format is not supported (code ",       &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_INDC_ADDR ) then
        print 33, "Error: array of indices is not defined (code ",       &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_WEIGHTS ) then
        print 33, "Error: array of weights contains negative values",    &
     &            "(code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MEAN_ADDR ) then
        print 33, "Error: array of means is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_2R_MOM_ADDR ) then
        print 33, "Error: array of 2nd order raw moments is not",        &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_3R_MOM_ADDR ) then
        print 33, "Error: array of 3rd order raw moments is not",        &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_4R_MOM_ADDR ) then
        print 33, "Error: array of 4th order raw moments is not",        &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_2C_MOM_ADDR ) then
        print 33, "Error: array of 2nd order central moments is not",    &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_3C_MOM_ADDR ) then
        print 33, "Error: array of 3rd order central moments is not",    &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_4C_MOM_ADDR ) then
        print 33, "Error: array of 4th order central moments is not",    &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_KURTOSIS_ADDR ) then
        print 33, "Error: array of kurtosis values is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_SKEWNESS_ADDR ) then
        print 33, "Error: array of skewness values is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MIN_ADDR ) then
        print 33, "Error: array of minimum values is not defined",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MAX_ADDR ) then
        print 33, "Error: array of maximum values is not defined",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_VARIATION_ADDR ) then
        print 33, "Error: array of variation coefficients is not",       &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_COV_ADDR ) then
        print 33, "Error: covariance matrix is not defined (code ",      &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_COR_ADDR ) then
        print 33, "Error: correlation matrix is not defined (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_QUANT_ORDER_ADDR ) then
        print 33, "Error: array of quantile orders is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_QUANT_ORDER ) then
        print 33, "Error: bad value of quantile order (code ",num,")."    
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_QUANT_ADDR ) then
        print 33, "Error: array of quantiles is not defined (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ORDER_STATS_ADDR ) then
        print 33, "Error: array of order statistics is not defined",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_MOMORDER_NOT_SUPPORTED ) then
        print 33, "Error: moment of requested order is not supported",   &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_NOT_FULL_RANK_MATRIX ) then
        print 33, "Warning: correlation matrix is not of full rank",     &
     &            " (code ",num,")."
      end if

      if ( num == VSL_SS_ERROR_ALL_OBSERVS_OUTLIERS ) then
        print 33, "Error: all observations are outliers (code ",         &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_COV_ADDR ) then
        print 33, "Error: robust covariance matrix is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_MEAN_ADDR ) then
        print 33, "Error: array of robust means is not defined",         &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_METHOD_NOT_SUPPORTED ) then
        print 33, "Error: requested method is not supported (code ",     &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_NULL_TASK_DESCRIPTOR ) then
        print 33, "Error: task descriptor is null (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OBSERV_ADDR ) then
        print 33, "Error: dataset matrix is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_SINGULAR_COV ) then
        print 33, "Error: covariance matrix is singular (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_POOLED_COV_ADDR ) then
        print 33, "Error: pooled covariance matrix is not defined",      &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_POOLED_MEAN_ADDR ) then
        print 33, "Error: array of pooled means is not defined (code ",  &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_COV_ADDR ) then
        print 33, "Error: group covariance matrix is not defined",       &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_MEAN_ADDR ) then
        print 33, "Error: array of group means is not defined (code ",   &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_INDC_ADDR ) then
        print 33, "Error: array of group indices is not defined (code ", &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_GROUP_INDC ) then
        print 33, "Error: group indices have improper values (code ",    &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_ADDR ) then
        print 33, "Error: array of parameters for outliers detection",   &
     &            " algorithm is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_N_ADDR ) then
        print 33, "Error: pointer to size of parameter array for",       &
     &            " outlier detection algorithm is not defined",         &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_OUTLIERS_WEIGHTS_ADDR ) then
        print 33, "Error: output of the outlier detection algorithm",    &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_ADDR ) then
        print 33, "Error: array of parameters of robust covariance",     &
     &            " estimation algorithm is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_N_ADDR ) then
        print 33, "Error: pointer to number of parameters of",           &
     &            " algorithm for robust covariance is not defined",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STORAGE_ADDR ) then
        print 33, "Error: pointer to variable that holds storage",       &
     &            " format is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COV_IDX_ADDR ) then
        print 33, "Error: array that encodes sub-components of",         & 
     &            " random vector for partial covariance algorithm",     &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COV_ADDR ) then
        print 33, "Error: partial covariance matrix is not defined",     &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COR_ADDR ) then
        print 33, "Error: partial correlation matrix is not defined",    &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PARAMS_ADDR ) then
        print 33, "Error: array of parameters for Multiple Imputation",  &
     &            " method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PARAMS_N_ADDR ) then
        print 33, "Error: pointer to number of parameters for",          &
     &            " Multiple Imputation method is not defined",          &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_BAD_PARAMS_N ) then
        print 33, "Error: bad size of the parameter array of Multiple",  &
     &            " Imputation method (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PARAMS ) then
        print 33, "Error: bad parameters of Multiple Imputation",        &
     &            " method (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_N_ADDR ) then
        print 33, "Error: pointer to number of initial estimates in",    &
     &            " Multiple Imputation method is not defined",          &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_ADDR ) then
        print 33, "Error: array of initial estimates for Multiple",      &
     &            " Imputation method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_SIMUL_VALS_ADDR ) then
        print 33, "Error: array of simulated missing values in",         &
     &            " Multiple Imputation method is not defined",          &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N_ADDR ) then
        print 33, "Error: pointer to size of the array of simulated",    &
     &            " missing values in Multiple Imputation method",       &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_ESTIMATES_N_ADDR ) then
        print 33, "Error: pointer to the number of parameter",           &
     &            " estimates in Multiple Imputation method is not",     &
     &            " defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_ESTIMATES_ADDR ) then
        print 33, "Error: array of parameter estimates in Multiple",     &
     &            " Imputation method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N ) then
        print 33, "Error: bad size of the array of simulated values",    &
     &            " in Multiple Imputation method (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_OUTPUT_PARAMS ) then
        print 33, "Error: array of output parameters in Multiple",       &
     &            " Imputation method is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PRIOR_N_ADDR ) then
        print 33, "Error: pointer to the number of prior parameters",    &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PRIOR_ADDR ) then
        print 33, "Error: bad number of missing values (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_MI_PRIOR_ADDR ) then
        print 33, "Error: bad number of missing values (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_SEMIDEFINITE_COR ) then
        print 33, "Warning: correlation matrix passed into",             &
     &            " parametrization function is semidefinite (code ",    &
     &             num,")."
      end if

      if ( num == VSL_SS_ERROR_BAD_PARAMTR_COR_ADDR ) then
        print 33, "Error: correlation matrix to be parametrized is",     &
     &            " not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_COR ) then
        print 33, "Error: all eigenvalues of correlation matrix to",     &
     &            " be parametrized are non-positive (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N_ADDR ) then
        print 33, "Error: pointer to the number of parameters for",      &
     &            " quantile computation algorithm for streaming data",  &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_ADDR ) then
        print 33, "Error: array of parameters of quantile computation",  &
     &            " algorithm for streaming data is not defined",        &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N ) then
        print 33, "Error: bad number of parameters of quantile",         &
     &            " computation algorithm for streaming data",           &
     &            " (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS ) then
        print 33, "Error: bad parameters of quantile computation",       &
     &            " algorithm for streaming data (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER_ADDR ) then
        print 33, "Error: array of quantile orders for streaming",       &
     &            " data is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER ) then
        print 33, "Error: bad quantile order for streaming data (code ", &
     &             num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_STREAM_QUANT_ADDR ) then
        print 33, "Error: array of quantiles for streaming data",        &
     &            " is not defined (code ",num,")."
        stop 1
      end if

      if ( num == VSL_SS_ERROR_BAD_PARTIAL_COV_IDX ) then
        print 33, "Error: partial covariance indices have improper",     &
     &            " values (code ",num,")."
        stop 1
      end if

33    format(A,I5,A)

      END SUBROUTINE
    
    end module eo_info
