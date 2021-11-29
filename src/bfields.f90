module bfields
      use cyclone_lib
      use cyclone_data
      implicit none      
      private

! public routines
      public BZIN,START,READ_CONE_FILE,READ_BUMP_FILE,BFLDCAL,bfldcal_odd,Harmonics,ModField,BZOUT,clean_bfields

! public variables
      real,   public,save  :: BCON,deltr,r0,rmx
      integer,public,save  :: nr,ntheta,n360
      REAL*8, public,save  :: DTHETA
      real,   public,save  :: BZ,dBZdT,dBZdR,d2BZdT2,d2BZdR2,BR,BTH
      real,   public,save  :: dBZdZ,d2BZdTdZ,d2BZdRdZ
! cone field
      real,   public,save  :: cone
! bump field
      real,   public,save  :: bump_mag,bump_ang
      integer,public,save  :: ibump_harm

! arrays local to this module
! main magnetic field
      real, save,allocatable :: ABZ(:,:)
      real, save,allocatable :: AdBdT(:,:)
      real, save,allocatable :: AdBdR(:,:)
      real, save,allocatable :: Ad2BdT2(:,:)
      real, save,allocatable :: ABR(:,:)
      real, save,allocatable :: ABTH(:,:)
      real, save,allocatable :: AdBdZ(:,:)
! cone field
      real, save,allocatable :: cone_field(:)
! bump field
      real, save,allocatable :: bump_field(:)

! generic interface for BFLDCAL
!  generic call is bfldcal(theta,r,reg_step,markz,*)
!  where theta may be either an integer or a real

      interface BFLDCAL
         module procedure integer_bfldcal
         module procedure real_bfldcal
      end interface

contains ! BZIN, read_cone_file,read_bump_file,START,FLDCAL,bfldcal_odd

! This subroutine creates (reads) the magnetic field arrays when the data
!  is stored in a file containing BZ values in an fixed r,theta grid
! The arrays store the field divided by BCON
! Use BFLDCAL to find field values at a point

      SUBROUTINE BZIN(name,format,file_header_lines)
! dummy variables
      CHARACTER*(*), Intent(IN) :: NAME   !File name for inout file
      CHARACTER*(*), Intent(IN) :: FORMAT ! Format string for reading the file - from input deck
      Integer, Intent(IN) :: file_header_lines ! true if the field file starts with a text header ended by EOF
      
! Local variables
 
      Character*csize actual_name,ErrMsg
      logical read_field_file
      logical file_header ! true if the field file starts with a text header ended by EOF
      
      Integer :: i,j,k,ik,ichar,NT,ir,ith,ierr
      Real :: R,theta,cos_theta

 300  FORMAT(A10)
 303  FORMAT(5I5,3F10.5)

!     In most cases field is read from file

      read_field_file=.true.
      file_header = .false.
      if (file_header_lines .gt. 0) file_header = .true.

! Allocate arrays

        if(allocated(abz)) deallocate(abz)
        allocate (ABZ(nr,ntheta+1),stat=ierr)
        call allocate_check(ierr,'Magnetic Field Array ABZ')        
        if(allocated(AdBdT)) deallocate(AdBdT)
        allocate (AdBdT(nr,ntheta+1),stat=ierr)
        call allocate_check(ierr,'Magnetic Field Array AdBdT ')        
        if(allocated(AdBdR)) deallocate(AdBdR)
        allocate (AdBdR(nr,ntheta+1),stat=ierr)
        call allocate_check(ierr,'Magnetic Field Array AdBdR')

!************* There is an error here in that we are testing against a fixed number of characters but name is dynamic in length
!
! Special case - B is uniform everywhere
      if(name(1:4) .eq. 'FLAT') then
      write(6,'("Using BFLAT")')
        do j=1,nr
        do i=1,ntheta
          abz(j,i)=BCON
        enddo
        enddo
        actual_name='FLAT'
        read_field_file=.false.
        goto 10
      endif

! Special Case - B is uniform in theta but has B=gamma*B0 shape in radius

      if(name(1:5) .eq. 'GAMMA') then
        ACON=E0*1000./(CHG*bconv*BCON*299.7925*conv)
        ASQ=ACON**2
        if (ASQ .eq. 0.0) then
          return_msg='ASQ cannot equal zero'
          ierror=1
          return
          !stop 'ASQ cannot equal zero'
        endif
        do j=1,nr
           r=r0+deltr*(j-1)
           do i=1,ntheta
              abz(j,i)=BCON/sqrt(1-r**2/asq)
           enddo
        enddo
        actual_name='GAMMA'
        read_field_file=.false.
        goto 10
      endif

10 continue
! General case where BZ is read from unit 44

      if(read_field_file) then
        OPEN (UNIT=44,STATUS='OLD',READONLY,file=name,err=99,iomsg=ErrMsg)
        inquire(44,name=actual_name)
        if(file_header) CALL START(file_header_lines,-44,0)

        DO J=1,NR
           READ (44,format,err=99,end=99,iomsg=ErrMsg),(ABZ(J,I),I=1,ntheta)
        end do
        close(44)
      endif
      
      ichar = len_trim(actual_name)
      write(6,'('' Magnetic Field from: '',a<ichar>)')actual_name(1:ichar)
      write(34,'(a<ichar>)')actual_name(1:ichar)
      write(6,'('' dr= '',f5.2,'' R0= '',f5.2,'' Nr= '',i5,"  B(0,0) =",f7.2)')deltr,r0,nr,abz(1,1)
      write(6,'(l5,i5)')file_header,file_header_lines

! Add Cone field (central field bump) if cone .ne. 0

      if(cone .ne. 0) then
        do j=1,nr
        do i=1,ntheta
          abz(j,i)=abz(j,i)+cone*cone_field(j)
        enddo
        enddo
      endif

! Add harmonic bump field if bump .ne. 0
!  for now it is assumed to be a first harmonic

      if(bump_mag .ne. 0) then
        if (nsc .ne. 1) then
           nt=nsc*ntheta
           IF(NT .GT. ntheta+1) THEN
              WRITE(6,'('' TOO MANY THETA VALUES '',2I5)')NT,ntheta+1
              return_msg='MAG. FIELD Bump LOADING ERROR'
              ierror=1
              return
              !STOP '*** MAG. FIELD Bump LOADING ERROR ***'
           END IF
           do ir=1,nr
              do ith=1,ntheta
                 do ik=1,nsc-1
                    abz(ir,ik*ntheta+ith)=abz(ir,ith)
                 enddo
              enddo
           enddo
           nsc=1
           ntheta=nt
        endif
        do i=1,ntheta
           theta=dtheta*(i-1)-bump_ang
           cos_theta=cos(ibump_harm*theta)
           do j=1,nr
              abz(j,i)=abz(j,i)+bump_mag*bump_field(j)*cos_theta
           enddo
        enddo
      endif

! DIVIDE FIELD BY BCON
!
      DO J=1,NR
      DO I=1,ntheta
        ABZ(j,i)=ABZ(j,i)/BCON
      end do
      end do
      call ComputeDerviatives
      RETURN 
      
 99   write(6,'('' Error reading field file '',A<csize>)',err=100)actual_name
      write(6,'('' Check format statement '',A20, ''  line= '',i5)',err=100)format,J
      write(6,'(A<csize>)',err=100)ErrMsg
      write(6,'("NR= ",i5," NTHETA= ",I5)',err=100)nr,ntheta
100   ierror=1
      return_msg=ErrMsg
      !Stop '*** Field Error***'
      return
      END SUBROUTINE BZIN
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Output BZ to a file
!
    SUBROUTINE BZOUT(name,format)

      CHARACTER*(*), Intent(IN) :: NAME   !File name for output file
      CHARACTER*(*), Intent(IN) :: FORMAT ! Format string for reading the file - from input deck

! local variables      
      real, allocatable :: Bout(:,:)
      Character*csize ErrMsg
      Integer :: i,j,ierr
      
! Allocate workspace
      allocate (Bout(nr,ntheta+1),stat=ierr)
! Stored field is divided by B0      
      Bout = ABZ*BCON

        OPEN (UNIT=44,STATUS='unknown',file=name)
        DO J=1,NR
           Write (44,format,err=99,iomsg=ErrMsg),(Bout(J,I),I=1,ntheta)
        end do
        close(44)
        write (6,'(/,"Wrote Isochronous field to ",A<csize> )')name
      
! Return allocated scratch space
      if(allocated(Bout)) deallocate(Bout)
      return
      
 99   write(6,'('' Error writing field file '',A<csize>)')name
      write(6,'('' Check format statement eg. (8F9.5) '')')
      write(6,'(A<csize>)')ErrMsg
      ierror=1
      return_msg=ErrMsg
      !Stop '*** Field Error***'
      return

      
    end subroutine BZOUT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Compute the dervitive of the main field
    Subroutine ComputeDerviatives

    ! Local variables
    integer :: i,j,k
!*** Should really have the computiuon of the derviatives in a single routine!!!!
    REAL*8 :: E1= 0.04166666667
    REAL*8 :: E2=-1.125d+0
    REAL*8 :: E3= 1.125d+0
    REAL*8 :: E4=-.04166666667

! dBZ/dTHETA
      do j=1,nr
        do i=2,ntheta-1
          AdBdT(j,i)=.5*(ABZ(j,i+1)-ABZ(j,i-1))/dtheta
        end do
        AdBdT(j,1)=.5*(ABZ(j,2)-ABZ(j,ntheta))/dtheta
        AdBdT(j,i)=.5*(ABZ(j,1)-ABZ(j,ntheta-1))/dtheta
      end do
! dBZ/dR
      K=NR-2 
      DO J=1,NTHETA 
        ADBDR(1,J)=(ABZ(2,J)-ABZ(1,J))/DELTR 
        ADBDR(K+1,J)=(ABZ(NR,J)-ABZ(NR-1,J))/DELTR 
        DO I=2,K 
          ADBDR(I,J)=(E1*ABZ(I-1,J)+E2*ABZ(I,J)+E3*ABZ(I+1,J)+E4*ABZ(I+2,J))/DELTR
        end do !i=2,k
      end do !j=1,ntheta
      return
      end subroutine ComputeDerviatives
             
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_cone_file(name,format,file_header)

      CHARACTER*(*) NAME
      CHARACTER*(*) FORMAT
      logical file_header 

      Character*csize actual_name
      Integer :: ichar,i,ierr

      OPEN (UNIT=45,STATUS='OLD',READONLY,file=name)      
      inquire(45,name=actual_name)
      ichar = len_trim(actual_name)
      if(file_header) CALL START(50,-44,0)
      write(6,'('' Cone Field from: '',a<ichar>)')actual_name(1:ichar)
      write(6,'('' Cone Magnitude = '',f8.5)')cone
      if(allocated(cone_field)) deallocate(cone_field)
      allocate (CONE_FIELD(nr),stat=ierr)
      call allocate_check(ierr,'Magnetic Field Array CONE_FIELD')        
      READ (45,format,err=99),(cone_field(i),i=1,nr)
      close(45)
      return
 99   write(6,'('' Cone file read error, nr='',i5)')nr
      write(6,'('' Format used: '',a10)')format
      return_msg='Cone file read error'
      ierror=1
      !stop
      return
      end subroutine read_cone_file
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_bump_file(name,format,file_header)

      CHARACTER*(*) NAME
      CHARACTER*(*) FORMAT
      logical file_header
      Integer :: Ierr,i,ichar

      Character*csize actual_name

      OPEN (UNIT=45,STATUS='OLD',READONLY,file=name)      
      inquire(45,name=actual_name)
      ichar = len_trim(actual_name)
      if(file_header) CALL START(50,-44,0)
      write(6,'('' Bump Field from: '',a<ichar>)')actual_name(1:ichar)
      write(6,'('' Bump Magnitude = '',f8.5)')bump_mag
      write(6,'('' Bump Angle = '',f8.3)')bump_ang*tcon
      write(6,'('' Bump Field harmonic = '',i3)')ibump_harm
      if(allocated(bump_field)) deallocate(bump_field)
      allocate (BUMP_FIELD(nr),stat=ierr)
      call allocate_check(ierr,'Magnetic Field Array BUMP_FIELD')        
      READ (45,format,err=99),(bump_field(i),i=1,nr)
      close(45)
      return
 99   write(6,'('' Bump file read error, nr='',i5)')nr
      write(6,'('' Format used: '',a10)')format
      return_msg='Bump file read error'
      ierror=1
      return
      !stop
      end subroutine read_bump_file
!**********************************************************************
      SUBROUTINE START(lines,IU,IOUT)
      
      integer, intent(in):: lines
      integer, intent(in):: IU
      integer, intent(in):: IOUT
            
! OPENS FILE ON UNIT IUNIT
! INQUIRES TO FIND THE FILE ON THAT UNIT
! READS UP TO N RECORDS AT BEGINNING OF FILE AND PRINTS THEM
! AS A DESCIPTION OF THE FILE

! N is the max number of lines to read from file
! IU IS THE INPUT FILE UNIT NUMBER
! IOUT IS THE UNIT NUMBER to WRITE THE INFO

! IF IU > 0 then open the file readonly, otherwise expect it to be already opened

! IF IOUT IS > 100 PRINTS INFO AND WRITES ON UNIT IOUT-100
! IF IOUT IS 0 DOES NOT WRITE INFO, JUST PRINTS INFO

! MODIFIED JUNE 10 1988
! MODIFIED sept 12 1988
! MODIFIED Jan 27 1999

	 LOGICAL*1 CHAR(80),IPRINT
     CHARACTER*csize FNAME
     Integer :: Iunit,ichar,i,j,nc,nrval,IO
         
100      FORMAT(Q,80A1)
101      FORMAT('  DATA FILE ',I5,'  IN ERROR')
102      FORMAT(2X,80A1)
103      FORMAT(' ',A<ichar>)
104      FORMAT(' ON UNIT ',I2,' IS FILE ',A<ichar>)

! open input file and obtain full file name

      IUNIT=ABS(IU)
      IF(IU .GT. 0) OPEN(IUNIT,STATUS='OLD',READONLY)
      INQUIRE(IUNIT,NAME=FNAME)

! remove trailing blanks from file name

      ichar=csize
        do while(fname(ichar:ichar) .eq. ' ')
        ichar=ichar-1
      end do

      IF(IOUT .GE. 100) THEN
        IPRINT=.TRUE.
        IO=IOUT-100
      else
        IPRINT=.FALSE.
        IO=IOUT
      END IF
      IF(IPRINT) PRINT 103,FNAME(1:ichar)
      IF(IO .NE. 0) WRITE(IO,104)IUNIT,FNAME(1:ichar)

! read in up to n lines from inuit and output as requested

      DO 10 I=1,abs(lines)
        READ(IUNIT,100,END=30,ERR=20),NC,(CHAR(J),J=1,80)
        IF(NC .GT. 80) NC=80
        IF(IPRINT)   PRINT 102,(CHAR(J),J=1,NC)
        IF(IO .NE.0) WRITE(IO,102),(CHAR(J),J=1,NC)
10    CONTINUE
      if (lines .gt. 0) return
! negative line count means EOF must be found

! error return
20    PRINT 101,IUNIT
      return_msg='FIELD FILE HEADER ERROR'
      ierror=1
      return
      !STOP ' *** FILE HEADER ERROR ***'

! normal exit

30    RETURN
      END SUBROUTINE START
!*********************************************************************** 
      LOGICAL FUNCTION REAL_BFLDCAL(theta,r,reg_step,markz)
      IMPLICIT NONE
! dummy variables
      REAL*8, INTENT(IN) :: theta
      REAL*8, INTENT(IN) :: r
      LOGICAL, INTENT(IN) :: reg_step
      INTEGER, INTENT(IN) :: markz
! local variables
      integer :: nt

      if(reg_step) then
         nt=nint(theta/dtheta)+1
      endif

      REAL_BFLDCAL = FLDCAL(theta,nt,r,reg_step,markz)
      return
      END FUNCTION REAL_BFLDCAL
!*********************************************************************** 
      LOGICAL FUNCTION INTEGER_BFLDCAL(ith,r,reg_step,markz)
      IMPLICIT NONE
! dummy variables
      INTEGER, INTENT(IN) :: ith
      REAL*8, INTENT(IN) :: r
      LOGICAL, INTENT(IN) :: reg_step
      INTEGER, INTENT(IN) :: markz
! local variables
      REAL*8 :: theta
      integer :: nt

      nt=ith
      theta=0
      if (.not. reg_step)then
         write(*,*) 'BFLDCAL called with integer theta, and reg_step=false'
         stop 'Illegal Call to BFLDCAL'
      endif

      INTEGER_BFLDCAL = FLDCAL(theta,nt,r,reg_step,markz)
      return

      END FUNCTION INTEGER_BFLDCAL
!*********************************************************************** 
! 
!  BFLDCAL:  COMPUTE FIELD AND DERIVATIVES
! 
!  theta = angle of field point in real*8
!  nt = angle specified as an integer
!  r = radius of field point in real*8
!  reg_step = logical switch - .true. if on a grid line
!  markz = order of calculation
!        = 0 - no z motion
!        = 1 - first order z motion
!        = 2 - second order z motion
!        = 3 - first order z with asymetric components
!        = 4 - second order z with asymetric components
!
!  A false return value means outside the radial limits of the field
!
!*********************************************************************** 
      LOGICAL FUNCTION FLDCAL(theta,nt,r,reg_step,markz)
!      LOGICAL FLDCAL
! dummy variables
      REAL*8,  INTENT(IN) :: theta
      INTEGER, INTENT(INOUT) :: nt
      REAL*8,  INTENT(IN) :: r
      LOGICAL, INTENT(IN) :: reg_step
      INTEGER, INTENT(IN) :: markz


! local variables
      LOGICAL ZEO
      REAL, DIMENSION(12) :: COF
      Real :: FRB,THETD,FTH,DBDZ1,DBDZ2,DBDZ3,DBDZ4,THETAD,PAREN
      Integer :: NRB,NRVAL,NTP,NTM1,NTP1,NTP2,INTV,NRVM1,NRBD
      Real :: DINT ! function
      
! internal function
      DINT(frb)=FLOAT(INT(frb))


! DETERMINE R POSITION

 10   continue
      NRVAL=nr
      FRB=(R-R0)/DELTR 
      NRB=INT(FRB) 
      FRB=FRB-FLOAT(NRB) 
      IF(R .LT. R0 .OR. R .GT. RMX) THEN
        if(reg_step) then
           thetad=(nt-1)*dtheta*tcon
        else
           thetad=theta*tcon
        endif
        WRITE(6,100) THETAD,R,rmx
  100   FORMAT('***OFF FIELD***  THETA=',E16.8,', RADIUS=',E16.8,' RMX = ',E16.8)
        FLDCAL = .false.
        RETURN
      END IF
      if(reg_step) then
        if(nt .le. 0)nt=nt+ntheta
        nt=mod(nt-1,ntheta)+1
        NTP=nt
        FTH=0.0
      else
        nt=ntheta+1
! do theta interp and stuff result in ntheta+1
        call bfldcal_odd(theta,r,markz,fth,ntp)
      end if
      ZEO=MARKZ .EQ. 2 .OR. MARKZ .EQ. 4
      IF(ZEO) THEN
        NTM1=NTP-1
        NTP1=NTP+1
        NTP2=NTP+2
        IF(NTM1.LT.1) NTM1=NTHETA
        IF(NTP2 .GT. NTHETA) NTP2=NTP2-NTHETA
        IF(NTP1 .GT. NTHETA) NTP1=NTP1-NTHETA
      END IF
 
! INTERPOLATION OF B AND THE THETA DERIVATIVES 
 
      BR=0.
      BTH=0.
      dBZdZ=0.
      INTV=1 ! First interval
!      IF(R .GE.  R0+DELTR) INTV=2 
!      IF(R .GE. RMX-DELTR) INTV=3
      if (nrb .ge. 1) INTV=2 ! Middle Intervals
      if (nrb .gt. nr-3) INTV=3 ! Final Interval

      select case(intv)

      case(1)
! FIRST INTERVAL 
   21 BZ=ABZ(1,NT)+FRB*(ABZ(2,NT)-ABZ(1,NT)) 
      IF(MARKZ .EQ. 0) GO TO 14
      dBZdT=AdBdT(1,NT)+FRB*(AdBdT(2,NT)-AdBdT(1,NT)) 
      IF(MARKZ .LT. 2) GO TO 14 
      IF(ZEO)d2BZdT2=Ad2BdT2(1,NT)+FRB*(Ad2BdT2(2,NT)-Ad2BdT2(1,NT)) 
      IF(MARKZ .LT. 3) GO TO 14
      BR=ABR(1,NT)+FRB*(ABR(2,NT)-ABR(1,NT)) 
      BTH=ABTH(1,NT)+FRB*(ABTH(2,NT)-ABTH(1,NT)) 
      dBZdZ=ADBDZ(1,NT)+FRB*(ADBDZ(2,NT)-ADBDZ(1,NT)) 
      IF(MARKZ .LT. 4) GO TO 14
      d2BZdRdZ=(ADBDZ(2,NT)-ADBDZ(1,NT))/DELTR 
      DBDZ1=ADBDZ(1,NTM1)+FRB*(ADBDZ(2,NTM1)-ADBDZ(1,NTM1))
      IF(REG_STEP) THEN
        DBDZ2=DBZDZ
      ELSE
        DBDZ2=ADBDZ(1,NTP)+FRB*(ADBDZ(2,NTP)-ADBDZ(1,NTP))
      END IF
      DBDZ3=ADBDZ(1,NTP1)+FRB*(ADBDZ(2,NTP1)-ADBDZ(1,NTP1))
      DBDZ4=ADBDZ(1,NTP2)+FRB*(ADBDZ(2,NTP2)-ADBDZ(1,NTP2))
      CALL INTCF(FTH,COF)
      d2BZdTdZ=(COF(5)*DBDZ1+COF(6)*DBDZ2+COF(7)*DBDZ3+COF(8)*DBDZ4)/DTHETA
      GO TO 14 
 
      case(3)
! LAST INTERVAL 
 
   23 NRVM1=NRVAL-1
      BZ=ABZ(NRVM1,NT)+FRB*(ABZ(NRVAL,NT)-ABZ(NRVM1,NT)) 
      IF(MARKZ .EQ. 0) GO TO 14
      IF(ZEO)DBZDT=AdBdT(NRVM1,NT)+FRB*(AdBdT(NRVAL,NT)-AdBdT(NRVM1,NT)) 
      IF(MARKZ .LT. 2) GO TO 14
      d2BZdT2=Ad2BdT2(NRVM1,NT)+FRB*(Ad2BdT2(NRVAL,NT)-Ad2BdT2(NRVM1,NT)) 
      IF(MARKZ .LT. 3) GO TO 14
      BR=ABR(NRVM1,NT)+FRB*(ABR(NRVAL,NT)-ABR(NRVM1,NT)) 
      BTH=ABTH(NRVM1,NT)+FRB*(ABTH(NRVAL,NT)-ABTH(NRVM1,NT)) 
      dBZdZ=ADBDZ(NRVM1,NT)+FRB*(ADBDZ(NRVAL,NT)-ADBDZ(NRVM1,NT))
      IF(MARKZ .LT. 4) GO TO 14
      d2BZdRdZ=(ADBDZ(NRVAL,NT)-ADBDZ(NRVM1,NT))/DELTR 
      DBDZ1=ADBDZ(NRVM1,NTM1)+FRB*(ADBDZ(NRVAL,NTM1)-ADBDZ(NRVM1,NTM1))
      IF(REG_STEP) THEN
        DBDZ2=DBZDZ
      ELSE
       DBDZ2=ADBDZ(NRVM1,NTP)+FRB*(ADBDZ(NRVAL,NTP)-ADBDZ(NRVM1,NTP))
      END IF
      DBDZ3=ADBDZ(NRVM1,NTP1)+FRB*(ADBDZ(NRVAL,NTP1)-ADBDZ(NRVM1,NTP1))
      DBDZ4=ADBDZ(NRVM1,NTP2)+FRB*(ADBDZ(NRVAL,NTP2)-ADBDZ(NRVM1,NTP2))
      CALL INTCF(FTH,COF)
      d2BZdTdZ=(COF(5)*DBDZ1+COF(6)*DBDZ2+COF(7)*DBDZ3 &
           +COF(8)*DBDZ4)/DTHETA
      GO TO 14 
 
      case(2)
! MIDDLE INTERVALS 
 
   22 CALL INTCF(FRB,COF)
      BZ=COF(1)*ABZ(NRB,NT)+COF(2)*ABZ(NRB+1,NT) &
           +COF(3)*ABZ(NRB+2,NT)+COF(4)*ABZ(NRB+3,NT)
      IF(MARKZ .EQ. 0) GO TO 14
      dBZdT=COF(1)*AdBdT(NRB,NT)+COF(2)*AdBdT(NRB+1,NT)  &
           +COF(3)*AdBdT(NRB+2,NT)+COF(4)*AdBdT(NRB+3,NT) 
      IF(MARKZ .LT. 2) GO TO 14
      IF(ZEO)d2BZdT2=COF(1)*Ad2BdT2(NRB,NT)+COF(2)*Ad2BdT2(NRB+1,NT) &
           +COF(3)*Ad2BdT2(NRB+2,NT)+COF(4)*Ad2BdT2(NRB+3,NT) 
      IF(MARKZ .LT. 3) GO TO 14
      BR=COF(1)*ABR(NRB,NT)+COF(2)*ABR(NRB+1,NT)+COF(3)*ABR(NRB+2,NT) &
           +COF(4)*ABR(NRB+3,NT)
      BTH=COF(1)*ABTH(NRB,NT)+COF(2)*ABTH(NRB+1,NT) &
           +COF(3)*ABTH(NRB+2,NT)+COF(4)*ABTH(NRB+3,NT)
      dBZdZ=COF(1)*ADBDZ(NRB,NT)+COF(2)*ADBDZ(NRB+1,NT) &
           +COF(3)*ADBDZ(NRB+2,NT)+COF(4)*ADBDZ(NRB+3,NT)
      IF(MARKZ .LT. 4) GO TO 14
      d2BZdRdZ=(COF(5)*ADBDZ(NRB,NT)+COF(6)*ADBDZ(NRB+1,NT) &
           +COF(7)*ADBDZ(NRB+2,NT)+COF(8)*ADBDZ(NRB+3,NT))/DELTR
      DBDZ1=COF(1)*ADBDZ(NRB,NTM1)+COF(2)*ADBDZ(NRB+1,NTM1) &
           +COF(3)*ADBDZ(NRB+2,NTM1)+COF(4)*ADBDZ(NRB+3,NTM1)
      IF(REG_STEP) THEN
        DBDZ2=DBZDZ
      ELSE
        DBDZ2=COF(1)*ADBDZ(NRB,NTP)+COF(2)*ADBDZ(NRB+1,NTP) &
             +COF(3)*ADBDZ(NRB+2,NTP)+COF(4)*ADBDZ(NRB+3,NTP)
      END IF
      DBDZ3=COF(1)*ADBDZ(NRB,NTP1)+COF(2)*ADBDZ(NRB+1,NTP1) &
           +COF(3)*ADBDZ(NRB+2,NTP1)+COF(4)*ADBDZ(NRB+3,NTP1)
      DBDZ4=COF(1)*ADBDZ(NRB,NTP2)+COF(2)*ADBDZ(NRB+1,NTP2) &
           +COF(3)*ADBDZ(NRB+2,NTP2)+COF(4)*ADBDZ(NRB+3,NTP2)
      CALL INTCF(FTH,COF)
      d2BZdTdZ=(COF(5)*DBDZ1+COF(6)*DBDZ2 &
           +COF(7)*DBDZ3+COF(8)*DBDZ4)/DTHETA
 
      end select

! RADIAL DERIVATIVES 
 
14    IF(MARKZ .EQ. 0) GO TO 40
      INTV=1 
      IF(R.GE.R0+DELTR/2.0) INTV=2 
      IF(R.GE.R0+3.0*DELTR/2.0) INTV=3 
      IF(R.GE.RMX-3.0*DELTR/2.0) INTV=4 
      IF(R.GE.RMX-DELTR/2.0) INTV=5 

      select case(intv)

! FIRST INTERVAL (LINEAR EXTRAPOLATION OF DBDR) 
 
      case(1)
   31 FRB=((R0+DELTR/2.0)-R)/DELTR 
      FRB=FRB-DINT(FRB) 
      dBZdR=ADBDR(1,NT)-FRB*(ADBDR(2,NT)-ADBDR(1,NT)) 
      IF(.NOT. ZEO) GO TO 40
      d2BZdR2=(ADBDR(2,NT)-ADBDR(1,NT))/DELTR 
      GO TO 40 
 
! SECOND INTERVAL (LINEAR INTERPOLATION OF DBDR) 
 
      case(2)
   32 FRB=(R-(R0+DELTR/2.0))/DELTR 
      FRB=FRB-DINT(FRB) 
      dBZdR=ADBDR(1,NT)+FRB*(ADBDR(2,NT)-ADBDR(1,NT)) 
      IF(.NOT. ZEO) GO TO 40
      d2BZdR2=(ADBDR(2,NT)-ADBDR(1,NT))/DELTR 
      GO TO 40 
 
! SECOND TO LAST INTERVAL 
 
      case(4)
   34 FRB=(R-(RMX-3.0*DELTR/2.0))/DELTR 
      FRB=FRB-DINT(FRB) 
      dBZdR=ADBDR(NRVAL-2,NT)+FRB*(ADBDR(NRVAL-1,NT)-ADBDR(NRVAL-2,NT)) 
      IF(.NOT. ZEO) GO TO 40
      d2BZdR2=(ADBDR(NRVAL-1,NT)-ADBDR(NRVAL-2,NT))/DELTR 
      GO TO 40 
 
! LAST INTERVAL 
 
      case(5)
   35 FRB=(R-(RMX-DELTR/2.0))/DELTR 
      FRB=FRB-DINT(FRB) 
      dBZdR=ADBDR(NRVAL-1,NT)+FRB*(ADBDR(NRVAL-1,NT)-ADBDR(NRVAL-2,NT)) 
      IF(.NOT. ZEO) GO TO 40
      d2BZdR2=(ADBDR(NRVAL-1,NT)-ADBDR(NRVAL-2,NT))/DELTR 
      GO TO 40 
 
! MIDDLE INTERVALS 
 
      case(3)
   33 FRB=(R-R0)/DELTR-.5
      NRBD=INT(FRB)
      FRB=FRB-DINT(FRB) 
      CALL INTCF(FRB,COF) 
      dBZdR=COF(1)*ADBDR(NRBD,NT)+COF(2)*ADBDR(NRBD+1,NT) &
           +COF(3)*ADBDR(NRBD+2,NT)+COF(4)*ADBDR(NRBD+3,NT) 
      IF(.NOT. ZEO) GO TO 40
      d2BZdR2=(COF(5)*ADBDR(NRBD,NT)+COF(6)*ADBDR(NRBD+1,NT) &
           +COF(7)*ADBDR(NRBD+2,NT)+COF(8)*ADBDR(NRBD+3,NT))/DELTR 

      end select

40    IF(ZEO) PAREN=DBZDR/R+D2BZDR2+D2BZDT2/R/R

    if (debug8) then
        write(6,'(1p4e12.5)')Bz,dBZdr,dBZdt
    endif
      FLDCAL = .true.
      RETURN 

      END FUNCTION FLDCAL
!*****************************************************
!
! Calculate the magnetic field for nonstandard steps
!  This routine does the theta interpolation and
!  then stores the result in ntheta+1, which is then 
!  used by BFLDCAL. BFLDCAL does the r interpolation.
!
!*****************************************************
      subroutine bfldcal_odd(theta,r,markz,fth,ith)
! dummy variables
      REAL*8,  INTENT(IN) :: theta
      REAL*8,  INTENT(IN) :: r
      INTEGER, INTENT(IN) :: markz
      Real, INTENT(INOUT) ::FTH
      Integer, INTENT(INOUT) :: ith
! local variables
      Real :: cof(8)
      Real :: FR,FT
      Integer :: IR,IR2,NT,KP,KR,KT,I,J

      fth=theta/dtheta
      ith=fth
      fth=fth-ith
      ith=mod(ith-1,ntheta)+1
      if(ith .eq. 0) ith=ntheta
      fr=(r-r0)/deltr
      ir=fr
      fr=fr-.5
      ir2=fr
      if(ir .lt. 1) ir=1
      if(ir2 .lt. 1) ir2=1
      if(ir+3 .gt. nr) ir=nr-3
      if(ir2+3 .gt. nr) ir2=nr-3
      nt=ntheta+1
      ft=fth
      call intcf(ft,cof)
      do i=1,4
        kr=i+ir-1
        kp=i+ir2-1
        abz(kr,nt)=0.
        if(markz .ge. 1) then
           adbdr(kp,nt)=0.
           adbdt(kr,nt)=0.
        endif
        if(markz .ge. 2) then
           ad2bdt2(kr,nt)=0.
        endif
        if(markz .ge. 3) then
           abr(kr,nt)=0.
           abth(kr,nt)=0.
           adbdz(kr,nt)=0.
        endif
        do j=1,4
          kt=ith+j-1
          if(kt .gt. ntheta) kt=kt-ntheta
          abz(kr,nt)=abz(kr,nt)+cof(j)*abz(kr,kt)
          if(markz .lt. 1) go to 50
          adbdr(kp,nt)=adbdr(kp,nt)+cof(j)*adbdr(kp,kt)
          adbdt(kr,nt)=adbdt(kr,nt)+cof(j)*adbdt(kr,kt)
          if(markz .lt. 2) go to 50
          ad2bdt2(kr,nt)=ad2bdt2(kr,nt)+cof(j)*ad2bdt2(kr,kt)
          if(markz .lt. 3) go to 50
          abr(kr,nt)=abr(kr,nt)+cof(j)*abr(kr,kt)
          abth(kr,nt)=abth(kr,nt)+cof(j)*abth(kr,kt)
          adbdz(kr,nt)=adbdz(kr,nt)+cof(j)*adbdz(kr,kt)
50        continue
        end do !j=1,4
      end do !i=1,4
      return
      end subroutine bfldcal_odd 
!*****************************************************
!
! Calculate the harmonics of the magnetic field
!  This is used bythe isochronizing routine
!
!*****************************************************
    subroutine harmonics(nhm,bav,g,h)
! nhm is the number of harmonics
! bav, h, & g and the returned harmonic information
    integer, intent(in) :: nhm
    real,intent(inout) :: bav(nr+15)
    real,intent(inout) :: g(nr,nhm)
    real,intent(inout) :: h(nr,nhm)

! this is local work space
     real*8, allocatable :: c(:,:)
     real*8, allocatable :: s(:,:)
     real*8, allocatable :: gt(:)
     real*8, allocatable :: ht(:)

! local variables     
     integer :: i,j,n,ierr
     real*8 :: b1,dth,th,thj

! allocate local storage    
     allocate (c(ntheta,nhm),stat=ierr)
     allocate (s(ntheta,nhm),stat=ierr)
     allocate (gt(nhm),stat=ierr)
     allocate (ht(nhm),stat=ierr)
          
     dth=TPI/ntheta
     do i=1,ntheta
        th=dth*(i-1)
        do j=1,nhm
           thj=j*th
           c(i,j)=cos(thj)
           s(i,j)=sin(thj)
        enddo
     enddo

     do n=1,nr
        ht = 0.d+0
        gt = 0.d+0
        bav(n)= 0.d+0
        do i=1,ntheta
           b1=2.D+0*abz(n,i)*bcon
           bav(n)=bav(n) + abz(n,i)
           do j=1,nhm
              ht(j)=ht(j)+b1*c(i,j)
              gt(j)=gt(j)+b1*s(i,j)
           enddo
        enddo
        do j=1,nhm
           h(n,j)=ht(j)/ntheta
           g(n,j)=gt(j)/ntheta
        enddo
        bav(n) = bcon*bav(n)/ntheta
     enddo

! return local storage
    deallocate (ht)
    deallocate (gt)
    deallocate (c)
    deallocate (s)
    return
    end subroutine harmonics

!**************************************************************************************
! Add an averge field component to the magnietic field at each radius, constant in theta    
!  Called from the isochronizing module
!**************************************************************************************
    subroutine ModField(Change,BCON_Initial)
    real, intent(in) :: change(nr) ! required change in average bz value
    Real, intent(in) :: BCON_Initial

    ! Local variables
    integer :: i,j,k
    Real :: Delta

! Apply correction to average field by adding to each point in theta 

    do i= 1,nr
        delta = change(i)
        do j=1,ntheta
           abz(i,j)=(abz(i,j)*bcon_initial + delta)/bcon
        end do
    end do
    
! Need to update orbit constants since BCON has changed
    call SetOrbitConstants
! Update field derviative - same as when field is read in
    call ComputeDerviatives

    return
    end subroutine ModField
!**********************************************************      
      subroutine clean_bfields
      
        if(allocated(ABZ)) deallocate(ABZ)
        if(allocated(AdBdT)) deallocate(AdBdT)
        if(allocated(AdBdR)) deallocate(AdBdR)
        if(allocated(Ad2BdT2)) deallocate(Ad2BdT2)
        if(allocated(ABR)) deallocate(ABR)
        if(allocated(ABTH)) deallocate(ABTH)
        if(allocated(AdBdZ)) deallocate(AdBdZ)
        if(allocated(cone_field)) deallocate(cone_field)
        if(allocated(bump_field)) deallocate(bump_field)
     
      end subroutine
         
end module bfields
