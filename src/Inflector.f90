module inflector

  use bAxial
  use eAxial
  use cyclone_lib
  use cyclone_data
  use io

  implicit none
  private

! public routines in this module
  public PartZero,ConfigPartZero,SetPartZeroPrint,SetCollimator
  
! public part 0 storage

! module level shared storage
    integer,private,save :: nstep=-4 ! number of degrees per step
    integer,private,save :: nprint=10 ! number of steps between prints
    integer,private,save :: max_step=1000 ! number of inflector steps (limit)
    real,private,save :: theta_match ! angle (in x-y) at which we leave the inflector calculations (procced to part 1)
    real,private,save :: r_match ! the radius at which we start checking theta_match
    real,private,save :: psign ! initial beam direction (+/-1), negative is down
    real*8,private,save :: zint ! initial vertical position
    real*4,private :: E_print ! Energy to print on unit 35
    real*8,private,save :: zCol ! vertical position of collimator
    real*8,private,save :: zColAbs ! vertical position of collimator (absolue value)
    real*8,private,save :: xCol,yCol ! width & breadth of collimator
    real*8,private,save :: thCol ! orientation of collimator axis - angle between x and xCol (rads)
    logical,private,save :: CheckCol ! do colimator check
    logical,private,save :: ColExists=.false. ! do colimator check    
  Contains

  SUBROUTINE PartZero(E,BCON,x_in,px_in,y_in,py_in,tau_in,theta_out)
    real*8, intent(inout) :: E ! energy
    real, intent(in) :: BCON ! central field
    real, intent(inout) :: x_in ! input x, output r
    real, intent(inout) :: px_in ! input px, output pr
    real, intent(inout) :: y_in ! input y, output z
    real, intent(inout) :: py_in ! input py, output pz
    real, intent(inout) :: tau_in ! input tau, output tau
    real, intent(out) :: theta_out ! return value for theta

!  we then figure out pz from the energy, z is assumed to be always the same (zint)  

    ! y(1) - tau
    ! y(2) - x
    ! y(3) - y
    ! y(4) - z
    ! y(5) - px
    ! y(6) - py
    ! y(7) - pz
    real*8 :: Y(7),Q(7)
    real*8 :: Efld(3) ! components of the electric field
    real*8 :: POT
    real*8 :: bx,by,bz
    real*8 :: PSQ,GAM,C1,C2,pc,gamma
    real*8 :: step
    integer :: istp ! number of integration steps - loop
    integer i,j,ierr
    Real :: Theta,Radius !end position

    PSQ=E/E0
    gamma=1.0D+0+psq
    psq=psq*(2.D+0+psq)*asq
    c1=chg*asq/e0
    c2=1.0/gamma
    pc=sqrt(psq)
    E_print=E
    checkCol=ColExists

! for positive nstep the step size will be nstep degs
! a negative nstep means that the step size will be 1/|nstep| degs
    if(nstep .ge. 0) then
      step=nstep*1.0D+0/tcon
    else
      step=-1.0D+0/tcon/nstep
    end if
! assume that if above MP then travel down, and if below then travel up    
    psign=sign(1.0,-zint)
           
    y(1)=tau_in/tcon !tau
    y(2)=x_in ! x
    y(3)=y_in ! y
    y(4)=zint  ! z
    y(5)=px_in ! px
    y(6)=py_in ! py
    y(7)=psign*sqrt(psq-px_in**2-py_in**2)
    do i=1,7
     q(i)=0.D+0
    end do

    istp=0

! NEED to figure out how to start above the electric field
!  this requires some concept of above (or below)

   
    call InputEaxial(psign)
    call InputBaxial(-Bcon) ! the minus sign is a hack that assumes the particle is negative
    
    call TitleEaxial ! print the electric field title information
    if(ColExists) write(6,'(" Inflector collimator at z=",f6.2,",",f6.2," by ",f6.2," at angle ",f7.2)')zcol,xcol,ycol,thcol*tcon
    write(6,102)

    call output(0,y)

    
! check that the ray is not already so large as to immediately transfer to part 1
    if(sqrt(y(2)**2+y(3)**2) .ge. r_match) then
      exitCode=iex_infstart
      write(6,'("Initial coordinates outside inflector x,y,r",3f8.3)')y(2),y(3),sqrt(y(2)**2+y(3)**2)
      go to 30
    endif
              
15  continue
      do j=1,4
        call getEaxial(Y,efld,pot,*80)
        CALL getBaxial(Y,BX,BY,BZ,*90)
        call intg(y,q,step,c1,c2,efld,bx,by,bz,j)
      end do
      !write(77,'(4f10.5,6e16.8)')y(1:4),efld,bx,by,bz
      istp=istp+1
      if(istp .gt. max_step) goto 20
      call output(istp,y)
      call check_match(Y,*20)
      if (CheckCol) call check_col(Y,*30)
    go to 15
    
20  continue
! now interpolate for theta_match

    radius=sqrt(y(2)**2+y(3)**2)
    theta=atan2(y(3),y(2))
    call thnrm(theta)
    !write(6,'("exit inflector theta,r,z",3f12.5)')theta*tcon,radius,y(4)
    write(6,'(" ")')
    call output(0,y)
    write(6,'(" ")')
        
    x_in=radius
    px_in=y(5)*cos(theta)+y(6)*sin(theta) ! pr
    y_in=y(4)  !z
    py_in=y(7) !pz
    theta_out=theta*tcon
    tau_in=y(1)*tcon

    return
    
30  continue ! hit collimator
        
80  continue ! off efield
    radius=sqrt(y(2)**2+y(3)**2)
    theta=atan2(y(3),y(2))
    call thnrm(theta)
    if (nlost .eq. 0) then
      ierr=openFile(85,return_msg)
    endif
    write(85,'(i5,5f12.5)')nrun,theta*tcon,radius,e,y(4),y(1)*tcon
    return    
    
90  continue !bfield error
    exitCode=iex_OffBfield
    return
    
! format statements
101   format(f8.3,3f8.4,3f10.5,f7.3,1x,1pe10.3,a1)
102   format(tr2'Tau',tr7,'X',tr7,'Y',tr7,'Z',tr8,'PX',tr8,'PY',&
       tr8,'PZ',tr6,'PT')         
  end subroutine
  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccc

! subroutine that does one step of a 4 part runge-kutta step

      subroutine intg(y,q,step,c1,c2,efld,bx,by,bz,j)
      implicit real*8 (a-h,o-z)
      real*8 y(7),ck(7),q(7),rrk,efld(3)
      real*8 ark(4)/.5D+0,.292893219D+0,1.707106781D+0,.1666666667D+0/
      real*8 brk(4)/2.D+0,1.D+0,1.D+0,2.D+0/
      real*8 crk(4)/-.5D+0,-.292893219D+0,-1.707106781D+0,-.5D+0/
      
      integer :: i,j

      ex=efld(1)
      ey=efld(2)
      ez=efld(3)
      ! find rrk coeficients
      ck(1)=1.0D+0
      ck(2)=c2*y(5)
      ck(3)=c2*y(6)
      ck(4)=c2*y(7)
      ck(5)=c1*EX+c2*(y(6)*bz-y(7)*by)
      ck(6)=c1*EY+c2*(y(7)*bx-y(5)*bz)
      ck(7)=c1*EZ+c2*(y(5)*by-y(6)*bx)
      do i=1,7
        ck(i)=ck(i)*step
      end do
      ! rrk loop
      do i=1,7
        rrk=ark(j)*(ck(i)-brk(j)*q(i))
        y(i)=y(i)+rrk
        q(i)=q(i)+3.D+0*rrk+crk(j)*ck(i)
      end do
      return
      end subroutine

!*****************************************
! Check if the particle has past the matching angle or not
!      
! need to carefully think about what the stop conditions are!
      subroutine check_match(Y,*)
          
          real*8,intent(in) :: Y(7)
          
          real theta ! angle in x-y plane
          real radius ! in the x-y plane
          real,save :: theta_previous
          
          radius=sqrt(y(2)**2+y(3)**2)
          if (radius .lt. r_match) then
            theta_previous=theta_match
            return
          endif
          
          theta=atan2(y(3),y(2))
          call thnrm(theta)
          if (theta .lt. theta_previous) theta_previous=theta_previous-TPI
          if (theta_previous .lt. 0.0 .and. theta .gt. 0.0) then ! straddles zero so should straddle 360 instead
            theta_previous=theta_previous+TPI
            theta=theta+TPI
          endif
          !write(6,'(4f12.5)')theta*tcon,theta_match*tcon,theta_previous*tcon,radius
          if (theta .ge. theta_match .and. theta_previous .lt. theta_match) return 1
          theta_previous=theta
          return
      
      end subroutine
!*****************************************
! Check if the particle has hit the collimator
!      
      subroutine check_col(Y,*)
          
          real*8,intent(inout) :: Y(7)
          
          real*8,save :: zprevious,xprevious,yprevious ! previous setup locations
          real*8 :: scale
          real*8 :: xc,yc ! values at collimator
          real :: xp,yp ! values in the collimator coordinate system

          if(abs(y(4)) .gt. zColAbs) then ! if before col then return
            xprevious=y(2)
            yprevious=y(3)
            zprevious=y(4)
            return
          endif
          
          ! check that not already past collimator
          if (abs(zprevious) .lt. abs(y(4))) then
            checkCol=.false.
            write(6,'(" It appears that this ray started after the collimator location so no collimator check",f8.4)')zprevious
            return
          endif
          
          ! now below collimator so interpolate for collimator
          scale=(zCol-zprevious)/(y(4)-zprevious)
          xc=xprevious+scale*(y(2)-xprevious)
          yc=yprevious+scale*(y(3)-yprevious)
          xp=xc*cos(thCol)+yc*sin(thCol)
          yp=-xc*sin(thCol)+yc*cos(thCol)
          write(10,'(i5,4f12.5)')nrun,xc,yc,xp,yp
          if((abs(xp) .ge. xCol) .or. (abs(yp) .ge. yCol)) then
            y(2)=xc
            y(3)=yc
            y(4)=zcol
            exitCode=iex_infcol
            write(6,'(" Particle hit collimator x,y ",2f8.4," Col Coord x,y",2f8.4)')xc,yc,xp,yp
            return 1
          endif
          checkCol=.false. ! stop checking
      end subroutine
      
!********************************************
    subroutine output(istp,y)
    real*8,intent(in) :: y(7)
    integer,intent(in) :: istp
    
    real*4 :: tau ! Y(1) in degrees
    real :: ptotal ! total momentum
    integer :: j ! general index
    real*4 :: theta,radius,pr ! values for output on unit 35 - precision must be maintained for binary file
    real*4 :: zero=0.0
    
      tau=y(1)*tcon
      ptotal=amp(y(5))

      if(mod(istp,nprint) .eq. 0) then
        write(6,101)tau,(y(j),j=2,7),ptotal
      endif
      IF(io_control35.active ) then
        io_control35.lines=io_control35.lines+1
        radius=sqrt(y(2)**2+y(3)**2)
        theta=atan2(y(3),y(2))
        call thnrm(theta)
        pr=y(5)*cos(theta)+y(6)*sin(theta) ! pr
         theta=theta*tcon
        if(io_control35.formatted) then
           WRITE(35,'(9e14.6,i2)')theta,radius,pr,0.0,e_print,TAU,y(4),y(7),0.,0
        else
           WRITE(35)theta,radius,pr,zero,e_print,TAU,real(y(4),4),real(y(7),4),zero,0.0
        endif
      endif ! active
      return

101   format(f8.3,3f8.4,3f10.5,f7.3,1x,1pe10.3,a1)   
    end subroutine output
    
!****************************************************
! Used by parameter load to get the PartZero settings
! called when Inflector command is found
!
    logical function ConfigPartZero
      real,external :: rdata_in
      
      zint=rdata_in(.false.,-1000.,1000.,0.)
      theta_match=rdata_in(.false.,0.,360.,0.)/tcon
      r_match=rdata_in(.false.,0.,1000.,1.)
   
      ConfigPartZero = .false.
    end function

! ***************************************************
!  called by parameter_load when print,0 is read
!    
    logical function SetPartZeroPrint
      integer,external :: idata_in
      
      nprint=idata_in(.false.,1,10000,1)
      nstep=idata_in(.false.,-100,1000,-4)
      max_step=idata_in(.false.,1,100000,1000)
      SetPartZeroPrint=.false. !no error
    end function SetPartZeroPrint

!****************************************************
! Used by parameter load to get the PartZero settings
! called when Collimator command is found
!
    logical function SetCollimator
      real,external :: rdata_in
      logical,external :: ldata_in
      
      ColExists=ldata_in(.false.,.true.)
      zCol=rdata_in(.false.,-1000.,1000.,0.)
      zColAbs=abs(zCol)
      xCol=rdata_in(.false.,0.,100.,1.)
      yCol=rdata_in(.false.,0.,100.,1.)
      thCol=rdata_in(.false.,0.,360.,1.)/tcon
      SetCollimator = .false.
      if(ColExists) then
        if (openFile(10,return_msg) .ne. 0)SetCollimator=.true.
        write(10,'(4f10.5)')zCol,xCol,yCol,thCol*tcon
      endif
    end function
end module