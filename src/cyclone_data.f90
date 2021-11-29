
! this module is the general common block for cyclone
!  used by most of the other cyclone modules

      module cyclone_data
      PUBLIC
      SAVE
      
! Error handling
      character*200 return_msg
      integer*4 ierror
      
! Physical Constants
      real*8, PARAMETER :: TPI=6.28318530717959D+0  !TWO PI
      real*8, PARAMETER :: TCON=360.0D+0/TPI        !CONVERT RADIANS TO DEGREES
      integer,parameter :: csize=200 ! this is the # of chars for filenames

! general orbit constants
      REAL*8  :: ASQ,VDEE,FRAT,AS,DCN1,ENPI,E0,CHG,ACON
      INTEGER :: NH,NTFIN,NSC,ndee
      Integer :: use_eps=0
      REAL    :: voltage,anurf,epsi,anu0,DIR
      REAL*8 :: comp_constant ! constant for momentum compaction

! Dee phasing and error specifciation
      REAL*8, allocatable :: DEEK(:),VERR(:),pherr(:),deekd(:)
      integer, allocatable :: iprn(:),idee(:)
      REAL*8, allocatable :: DEEK3(:),VERR3(:),pherr3(:)
      Logical :: third_harmonic_rf

! output control
      integer :: nprnt2,nprntd2,nprntt2,nprnt3,nprntd3,nprntt3,nrf
      logical :: lsect2,lsect3 ! repeat by sector
      integer :: nsect2,nsect3 ! number of steps in sect
      logical :: debug1,debug2,debug3,debug4,debug5,debug6,debug7,debug8

! counter of the number of particles run
     integer,public,save :: nrun !number of particles started
     
! Exit codes
    integer,parameter :: iex_nostrip=1
    integer,parameter :: iex_stripped=2
    integer,parameter :: iex_turncount=3
    integer,parameter :: iex_stoppt3=4
    integer,parameter :: iex_post=5
    integer,parameter :: iex_OffBfield=6
    integer,parameter :: iex_OffEfield=7
    integer,parameter :: iex_stoppt2=8
    integer,parameter :: iex_FoilFrame=9
    integer,parameter :: iex_CalcFailure=10
    integer,parameter :: iex_Undefined=11
    integer,parameter :: iex_WriteError=12
    integer,parameter :: iex_infpost_p=13
    integer,parameter :: iex_infcol=14 
    integer,parameter :: iex_infpost_n=15
    integer,parameter :: iex_infstart=16
    integer,parameter :: iex_FlagHit=17    
    integer :: ExitCode=0
    integer :: nlost=0 ! lost particle counter
! units
      character*2 :: units,bunits
      real :: conv,bconv
! transfer controls
      integer :: nturn_pt2,nturn_pt3
      real    :: theta_pt2,theta_pt3
      logical ::  stop_pt0,stop_pt2,stop_pt3

! input data control information
      logical :: ReadFile ! Either read from unit 5 or get information from InputLines
      Integer :: numInputLines ! number of rows in InputLines
      Character*csize, allocatable :: InputLines(:)
      Integer :: currentLine
      
      contains
      
      subroutine clean_cyclone_data
      
        if(allocated(deek)) deallocate(deek)
        if(allocated(VERR)) deallocate(VERR)
        if(allocated(pherr)) deallocate(pherr)
        if(allocated(deekd)) deallocate(deekd)
        if(allocated(iprn)) deallocate(iprn)
        if(allocated(idee)) deallocate(idee)
        if(allocated(deek3)) deallocate(deek3)
        if(allocated(VERR3)) deallocate(VERR3)
        if(allocated(pherr3)) deallocate(pherr3)
        if(allocated(InputLines)) deallocate(InputLines)
      
      end subroutine
      
      end module cyclone_data
