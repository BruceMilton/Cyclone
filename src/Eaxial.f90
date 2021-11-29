module Eaxial

    use cyclone_lib
    !use cyclone_data, only: ierror,return_msg
    use cyclone_data
    use IO
    implicit none

    private

! public routines in this module
    public SetE_Axial,InputEaxial,GetEaxial,TitleEaxial,Set_Axial_Posts

! public inflector (part 0) storage
    real,public,save :: InfVoltage
        
! module level shared storage

    real*4,private,save,allocatable :: potential(:,:,:)

    ! th_rot = ROTATION ANGLE OF GRID'S +X AXIS FROM THE THETA=0 LINE OF THE 
    !             MAGNETIC FIELD (COUNTERCLOCKWISE IS POSITIVE DIRECTION)
    real,private,save :: th_rot,thr_rot,cos_rot,sin_rot ! rotation angle from input deck
    logical,private,save :: rotate ! flag indicates field rotation is active

    logical,private,save :: eLoaded=.false. ! flag to indicate the eAxial is loaded
    real,private,save :: scale=1.0 ! linear scaling of grid dimensions
    real,private,save :: zmax,zmin ! upper and lower bounds of the efield
    logical,private,save :: checkAbove ! switch to all above or below the inflector to run
    
    real,private,save :: dxl,dyl,dzl ! grid spacing
    ! XP0 = DISTANCE FROM CYCLOTRON CENTER TO RIGHT SIDE OF GRID
    ! YP0 = DISTANCE FROM CYCLOTRON CENTER TO BOTTOM OF GRID
    ! ZP0 = Distance from k=1 plane to median plane
    real,private,save :: xp0,yp0,zp0
    ! NYL = NUMBER OF GRID POINTS IN Y DIRECTION
    ! NXL = NUMBER OF GRID POINTS IN X DIRECTION
    ! NZL = NUMBER OF Z PLANES
    integer,private,save :: nxl,nyl,nzl ! number of points
    
    integer,private,save :: max_hits ! maximum post count in an e field calculation
    logical,private,save :: no_posts_allowed=.false. ! enable post checking
    logical,private,save :: post_log=.false. ! print all post points on unit 6
       
    character*csize filename
    
contains

! Input the electric field for the inflector section
! The electric field is stored in one RELAX3D mesh 
   SUBROUTINE InputEaxial(psign)
      real*4,intent(in) :: psign ! direction of initial travel

      Character*csize file_name
      Character*csize fname
      character*csize errmsg ! error message from IO errors
      integer :: ijunk(5)
      integer :: iunit1,iunit3 ! unit numbers for head and efld files
      
      real :: xpi,ypi,zpi ! header file offsets
      real :: dxi,dyi,dzi ! header grid spacing
      integer :: nxi,nyi,nzi ! header grid sizes
      
      integer :: ipp,ip,ierr,i,j,k
      
      if(psign .lt. 0) then
        checkAbove=.true.
      else
        checkAbove=.false.
      endif

      if(eLoaded) return ! field is already loaded
      
      iunit1=98
      iunit3=99
      file_name=filename(1:len_trim(filename))//'.head'
      OPEN(iunit1,STATUS='OLD',READONLY,FORM='FORMATTED',FILE=file_name,err=98,iomsg=ErrMsg)
      file_name=filename(1:len_trim(filename))//'.efld'
      OPEN(iunit3,STATUS='OLD',READONLY,FORM='UNFORMATTED',ACCESS='Stream',FILE=file_name,err=98,iomsg=ErrMsg)
      inquire(iunit3,name=fname)
      write(6,'(''Axial Electric Field from file '',a132)')fname

      ! read grid dimensions from the header file

      read(iunit1,*)xpi,ypi,zpi
      read(iunit1,*)dxi,dyi,dzi
      read(iunit1,*)nxi,nyi,nzi

      ! now read the problem dimensions from the relax file and compare
      !  if different stop program

      ipp=(nxi+1)*(nyi+1)*(nzi+1)
      read(iunit3)nxl,nyl,nzl,ijunk
      ip=nxl*nyl*nzl
      if(ip .ne. ipp) then
        write(6,'('' *** The number of Points in data file does not agree with header data '')')
        write(6,'(''File '',4i6)')ip,nxl,nyl,nzl
        write(6,'(''Header '',4i6)')ipp,nxi+1,nyi+1,nzi+1
        ierror=20
        return_msg='Eaxial error stop'
        return
        !stop '*** Error stop ***'
      end if

      !  get space to store array

      if(allocated(potential)) deallocate(potential)
      allocate (potential(nxl,nyl,nzl),stat=ierr)
      call allocate_check(ierr,'Inflector Electric Field Array')      
      !write(6,'( '' Points '',3i5,i7)')nxl,nyl,nzl,ip

      !  input relax data

      read(iunit3,end=97,err=97,iomsg=ErrMsg) (((potential(i,j,k),i=1,nxl),j=1,nyl),k=1,nzl)

      ! close files and return the units

      CLOSE(iunit1)
      CLOSE(iunit3)

      ! it is assumed that the relaxation calc was done with one plate
      !  at -3 and the other at -1 potential. The ground planes would then
      !  be at -2. We now rescale the problem so that pot lies between -1 and +1.
      !potential = abs(potential)-2.0
      
      eLoaded=.true.

      IF(SCALE .EQ. 0.) SCALE=1.
      dxl=dxi*scale
      dyl=dyi*scale
      dzl=dzi*scale
      XP0=XPi*SCALE
      YP0=YPi*SCALE
      ZP0=ZPi*SCALE
      WRITE(6,4020),NXL,NYL,NZL,DXL,DYL,DZL,XP0,YP0,ZP0
4020  FORMAT(' Inflector Electric Field',' NX=',I4,' NY=',I4,' NZ=',I4,/,1X,' DX=',F7.4,' DY=',F7.4,' DZ=',F7.4,' XP0=',F7.3,' YP0=',F7.3,' ZP0=',F7.3)
      zmax=max(zp0,zp0+dzl*(nzl-1))
      zmin=min(zp0,zp0+dzl*(nzl-1))

      RETURN

! ERROR EXITS

98    write(6,'(" Error opening electric field file")')
      write(6,'(a200)')errmsg
      write(6,'(a200)')file_name
      stop '*** Error opening relax file ***'
      
97    write(6,'('' Problem with RELAX3D file '')')
      write(6,4030)nxl,nyl,nzl
4030  format(5x,'requested dimensions for relax',/,'nx = ',i5,/,'ny = ',i5,/,'nz = ',i5)
      STOP '*** relax file error ***'
      
      end subroutine InputEaxial

!**************************************************************

!  print the electric field title information on unit 6

  subroutine TitleEaxial
  
      ! need voltage in mega volts
      write(6,300)InfVoltage*1000.,th_rot
300   format(" The Inflector Plates are at +/- ",f7.3," kV,"," The grid is rotated by ",f8.3," degs")
  end subroutine
  
!**************************************************************  
! Calculate the electric field at th,r,z and rf time tau

! The electric field is stored in cartesian coordinates.
! The x axis lies along the theta=0 axis unless th_rot is used to shift it. 
!   th_rot is a counterclockwise rotation about the z axis. 
! The number of points in each dimension is given by NXl,NYl,NZl.
! The spacings are given by DXL,DYL,DZL.
! The origin of the grid relative to r=0 is given by xp0, yp0, and zp0
!  which are negative if r=0 lies inside the grid, ie. they are
!  the coordinates of (1,1,1) in a system centered at r=0.
! Double three point lagrange interpolation is used to calculate
!  the potential and the fields at the point with x=pos(1), y=pos(2),
!  and z=pos(3). At the edges of the grid a single three point interpolation is used.

    Subroutine GetEaxial(pos,efld,pot,*)
      real*8,intent(in) :: pos(*) ! tau,x,y,z location to find fields at
      real*8,intent(out) :: efld(3) ! electric field components
      real*8, intent(out) :: pot ! electric potential

      real :: COFX(12),COFY(12),COFZ(12)
      real :: x,y,z,xp,yp
      real :: fx,fy,fz ! fractional step
      integer :: ix,iy,iz ! grid index
      real :: d ! potential temp storage
      real*8 :: efx,efy,efz,dfx,dfy,dfz,ff,test,ex,ey,ez,ddx,ddy,ddz ! variables to collect the fields and derivatives
      integer :: imax,jmax,kmax ! loop limits when closer than 2 point to edge of grid
      integer :: i,j,k,lerror
      integer :: icount

      if (rotate) then
        x=cos_rot*pos(2)+sin_rot*pos(3)
        y=-sin_rot*pos(2)+cos_rot*pos(3)
        xp=x-xp0
        yp=y-yp0
      else
        xp=pos(2)-xp0
        yp=pos(3)-yp0
      endif
      z=pos(4)-zp0
      if(IO_Control53.active)write(53,'(3(1p,e16.8))')z,pos(4),zp0
      FX=XP/dxl
      FY=YP/dyl
      FZ=Z/DZL
      if(IO_Control53.active)write(53,'('' POS '',9f10.5)')pos(2:4),xp,yp,z,fx,fy,fz
      IX=FX
      IY=FY
      IZ=FZ
      FX=FX-IX
      FY=FY-IY
      FZ=FZ-IZ
      if(io_control53.active)write(53,'(3(e16.8,i5))')fx,ix,fy,iy,fz,iz
      IF(IX .LT. 1 .OR. IX+3 .GT. NXL) GO TO 190
      IF(IY .LT. 1 .OR. IY+3 .GT. NYL) GO TO 190
      IF(IZ .LT. 1 .OR. IZ+3 .GT. NZL) GO TO 190
      CALL INTCF12(FX,COFX)
      CALL INTCF12(FY,COFY)
      CALL INTCF12(FZ,COFZ)
      EFX=0.D+0
      EFY=0.D+0
      EFZ=0.D+0
      FF =0.D+0
      DDX=0.D+0
      DDY=0.D+0
      DDZ=0.D+0
      icount=0
      DO 50 K=1,4
      DO 50 J=1,4
      DO 50 I=1,4
        d=potential(ix+i-1,iy+j-1,iz+k-1)
        if(d .lt. 0) icount=icount+1
        d=abs(d)-2.0
        EFX=EFX-COFX(I+4)*COFY(J)*COFZ(K)*D
        EFY=EFY-COFY(J+4)*COFX(I)*COFZ(K)*D
        EFZ=EFZ-COFZ(K+4)*COFX(I)*COFY(J)*D
        FF=FF+COFX(I)*COFY(J)*COFZ(K)*D
        DDX=DDX+cofx(i+8)*cofy(j)*cofz(k)*d
        DDY=DDY+cofy(j+8)*cofx(i)*cofz(k)*d
        DDZ=DDZ+cofz(k+8)*cofx(i)*cofy(j)*d
50    CONTINUE !J,I,K
      efx=efx/dxl
      efy=efy/dyl
      ez =efz/dzl
      if(rotate) then
        EX=EFX*cos_rot-EFY*sin_rot
        EY=EFY*cos_rot+EFX*sin_rot
      else
        EX=EFX
        EY=EFY
      endif
      TEST=DDX/dxl/dxl+DDY/dyl/dyl+DDZ/dzl/dzl
      efld(1)=InfVoltage*EX
      efld(2)=InfVoltage*EY
      efld(3)=InfVoltage*ez
      POT=FF*InfVoltage
      if(io_control53.active) then
        write(53,'('' X '',4f12.7)')(potential(ix+i,iy+2,iz+2),i=0,3)
        write(53,'('' Y '',4f12.7)')(potential(ix+2,iy+i,iz+2),i=0,3)
        write(53,'('' Z '',4f12.7)')(potential(ix+2,iy+2,iz+i),i=0,3)
        write(53,'('' CX'',1p,8e16.8)')(cofx(i),i=1,8)
        write(53,'('' CY'',1p,8e16.8)')(cofy(i),i=1,8)
        write(53,'('' CZ'',1p,8e16.8)')(cofz(i),i=1,8)
        write(53,'(1p,6e16.8)')efx,efy,efz,ff,InfVoltage,test
        write(53,'('' EFLD,POT '',1p,4e16.8)')efld,pot
      end if
      goto 90

!  less than 2 points from edge

190   continue
      IF(FX .LT. 0.)GOTO 191
      IF(FY .LT. 0.)GOTO 191
      IF(FZ .LT. 0.)GOTO 191
      IF(IX .LT. 0 .OR. IX+2 .GT. NXL) GO TO 191
      IF(IY .LT. 0 .OR. IY+2 .GT. NYL) GO TO 191
      IF(Iz .LT. 0 .OR. Iz+2 .GT. NzL) GO TO 191
      if(ix .eq. 0) then
        fx=fx-1.
        CALL INTCF3D(FX,COFX)
        imax=3
        ix=ix+1
      else if(ix+2 .eq. nxl) then
        CALL INTCF3D(FX,COFX)
        imax=3
      else
        CALL INTCFD(FX,COFX)
        imax=4
      end if

      if(iy .eq. 0) then
        fy=fy-1.
        CALL INTCF3D(FY,COFY)
        jmax=3
        iy=iy+1
      else if(iy+2 .eq. nyl) then
        CALL INTCF3D(FY,COFY)
        jmax=3
      else
        CALL INTCFD(FY,COFY)
        jmax=4
      end if

      if(iz .eq. 0) then
        fz=fz-1.
        CALL INTCF3D(FZ,COFZ)
        kmax=3
        iz=iz+1
      else if(iz+2 .eq. nzl) then
        CALL INTCF3D(FZ,COFZ)
        kmax=3
      else
        CALL INTCFD(FZ,COFZ)
        kmax=4
      end if
      EFX=0.D+0
      EFY=0.D+0
      EFZ=0.D+0
      FF =0.D+0
      icount=0
      DO 70 K=1,kmax
      DO 70 J=1,jmax
      DO 70 I=1,imax
        d=potential(ix+i-1,iy+j-1,iz+k-1)
        if(d .lt. 0) icount=icount+1
        d=abs(d)-2.0
        EFX=EFX-COFX(I+4)*COFY(J)*COFZ(K)*D
        EFY=EFY-COFY(J+4)*COFX(I)*COFZ(K)*D
        EFZ=EFZ-COFZ(K+4)*COFX(I)*COFY(J)*D
        FF=FF+COFX(I)*COFY(J)*COFZ(K)*D
70    CONTINUE !J,I,K
      efx=efx/dxl
      efy=efy/dyl
      ez =efz/dzl
      if(rotate) then
        EX=EFX*cos_rot+EFY*sin_rot
        EY=EFY*cos_rot-EFX*sin_rot
      else
        EX=EFX
        EY=EFY
      endif
      efld(1)=InfVoltage*EX
      efld(2)=InfVoltage*EY
      efld(3)=InfVoltage*EZ
      POT=FF*InfVoltage
      if(io_control53.active) then
        write(53,'('' X '',4f12.7)')(potential(ix+i,iy+2,iz+2),i=0,imax-1)
        write(53,'('' Y '',4f12.7)')(potential(ix+2,iy+i,iz+2),i=0,jmax-1)
        write(53,'('' Z '',4f12.7)')(potential(ix+2,iy+2,iz+i),i=0,kmax-1)
        write(53,'('' CX'',8f10.5)')(cofx(i),i=1,8)
        write(53,'('' CY'',8f10.5)')(cofy(i),i=1,8)
        write(53,'('' CZ'',8f10.5)')(cofz(i),i=1,8)
        write(53,'(6f10.5)')efx,efy,efz,ff,InfVoltage,test
        write(53,'('' EFLD,POT '',1p,4e16.8,i5)')efld,pot,icount
      end if
90    continue ! back to common code for 4 and 3 point interp
      if(io_control62.active)write(62,'(5f10.5)')efld,amp(efld),test      
      lerror=0
      if(icount .gt. 0) Then
        if(post_log) write(6,'("Electrode IX,IY,IZ,V,Count",3i6,e12.3,i3)')ix,iy,iz,pot,icount
        if(no_posts_allowed .and. icount .ge. max_hits) then
          write(6,'("Particle lost due to number of electrode points used, X,Y,Z,Hits ",3f10.5,i5)')pos(2:4),icount
          if(pot .ge. 0) then
            exitCode=iex_infpost_p
          else
            exitCode=iex_infpost_n
          endif
          return 1
        endif
      endif
      return

! OFF FIELD

191   continue
      if(io_control53.active)write(53,'("Off Field ",3i5)')ix,iy,iz
      if(lerror.ne. 1) then
        WRITE(6,'('' OFF AXIAL E FIELD '',3i5)')ix,iy,iz
        write(6,'('' x = '',f7.3,'' y = '',f6.2,'' Z = '',f7.3)')xp,yp,z
        lerror=1
      end if
      efld(1)=0.0D+0
      efld(2)=0.0D+0
      efld(3)=0.0D+0
      pot=0.0D+0
      if(checkAbove) then
        if (pos(4) .gt. zmax) return ! above inflector & coming down
      else
        if (pos(4) .lt. zmin) return ! below inflector & coming up
      endif
      exitCode=iex_OffEfield
      return 1
    end Subroutine GetEaxial
    
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Called by Parameter_load to setup axial field read

  logical function SetE_Axial
    real*4,external :: rdata_in
  
    infVoltage=rdata_in(.false.,real(-1.0e6,4),real(1.0e6,4),real(0.0,4))/1000.
    th_rot=rdata_in(.false.,real(-360.0,4),real(360.0,4),real(0.0,4))
    if(th_rot .ne. 0.0) then
      cos_rot=cosd(th_rot)
      sin_rot=sind(th_rot)
      thr_rot=th_rot/tcon
      rotate=.true.
    else
      thr_rot=0.0
      cos_rot=1.0
      sin_rot=0.0
      rotate=.false.
    endif
    call getLine(filename,*10,*10)
    !read(5,'(a<csize>)')filename ! get field file name
    filename=adjustl(filename)
    eLoaded=.false.
    SetE_Axial=.false.
    return
10  SetE_Axial=.true.    
  end function SetE_Axial
  
! Called by Parameter_load to setup axial post checking

  logical function Set_Axial_Posts
    logical,external :: ldata_in
    integer,external :: idata_in
    
    max_hits=idata_in(.false.,0,100,4)
    no_posts_allowed=ldata_in(.false.,.true.)
    post_log=ldata_in(.false.,.false.)
    
    Set_Axial_Posts=.false.    
  end function Set_Axial_Posts
end module