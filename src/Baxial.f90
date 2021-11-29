module Baxial

    use cyclone_lib
    use cyclone_data, only: ierror,return_msg
    use IO    
    implicit none

    private

! public routines in this module
    public SetBaxial,InputBaxial,GetBaxial

! public part 1 storage

! module level shared storage
    real*4,private,save :: B0 ! save the value of BCON
    real*4,private,save,allocatable :: bx(:,:,:),by(:,:,:),bz(:,:,:),br(:,:,:),bth(:,:,:)
    integer,private,save :: nz,nr,nth,nx,ny !number of points
    integer,private,save :: ndim ! number of active dimensions
    real,private,save :: x0,y0,z0,r0,th0 ! initial grid values
    real,private,save :: dx,dy,dz,dr,dth ! grid spacing
    logical,private,save :: symx,symy,symz ! flags to indicate if these are symmetry planes
    real,private,save :: sym_ang  ! rotation symmetry about z for polar files - repeat angle in rads
    character*2,private,save :: dType
    real,private,save :: th_rot,thr_rot,cos_rot,sin_rot ! rotation angle from input deck
    logical,private,save :: rotate ! flag indicates field rotation is active
    logical,private,save :: bLoaded=.false. ! flag to indcate the bAxial is loaded
    
    integer,parameter :: csize=200 ! this is to become the #of chars for filenames
    character*csize filename
    
    real(8),parameter :: pi=3.141592653589793D+0
    real(8),parameter :: tcon=180.D+0/pi
    real(8),PARAMETER :: TPI=6.28318530717959D+0 !TWO PI
    real(8),parameter :: zero=0.0D+0
    real(8),parameter :: one=1.0D+0

  contains

  subroutine InputBaxial(BCON)
    real, intent(in) :: BCON ! central field
    
    character*5 ftype
    character*500 comment
    character*csize actual_name
    character*csize errmsg
    integer :: ierr
    real*4 :: x1,x2,x3,d1,d2,d3
    integer :: n1,n2,n3
    integer :: i
    
    if(bLoaded) return ! field is already input
    
    B0=BCON
    
    OPEN (UNIT=60,STATUS='OLD',READONLY,FORM='UNFORMATTED',file=filename,err=99,iomsg=errmsg)
    inquire(60,name=actual_name)
 
    read(60,err=97,iomsg=errmsg)ftype,ndim,n1,n2,n3
    read(60,err=98,iomsg=errmsg)x1,x2,x3,d1,d2,d3
    read(60,err=98,iomsg=errmsg)symx,symy,symz,sym_ang
    read(60,err=98,iomsg=errmsg)comment

    if(allocated(bz)) deallocate(bz)
    if(allocated(bx)) deallocate(bx)
    if(allocated(by)) deallocate(by)
    if(allocated(br)) deallocate(br)
    if(allocated(bth)) deallocate(bth)
    
    select case(ftype)
      Case ('POLAR')
        nz=n1
        nr=n2
        nth=n3
        z0=x1
        r0=x2
        th0=x3/tcon
        dz=d1
        dr=d2
        dth=d3/tcon
        write(6,'("Axial magnetic field from file ",a132)')actual_name
        write(6,'("z",i5,2f8.2)')nz,z0,dz
        write(6,'("r",i5,2f8.2)')nr,r0,dr
        write(6,'("th",i5,2f8.2)')nth,th0*tcon,dth*tcon
        select case(ndim)
          case(1)
            allocate (Bz(n1,1,1),stat=ierr)
            call allocate_check(ierr,'Bz polar')
            read(60,err=98,iomsg=errmsg)Bz
            Bz=Bz/bcon
            dtype='P1'
           case(2)
            allocate (Bz(n1,n2,1),stat=ierr)
            call allocate_check(ierr,'Bz polar')
            allocate (Br(n1,n2,1),stat=ierr)
            call allocate_check(ierr,'Br polar')
            read(60,err=98,iomsg=errmsg)Bz
            read(60,err=98,iomsg=errmsg)Br
            Bz=Bz/bcon
            Br=Br/bcon
            dtype='P2'
           case(3)
            allocate (Bz(n1,n2,n3),stat=ierr)
            call allocate_check(ierr,'Bz polar')
            allocate (Br(n1,n2,n3),stat=ierr)
            call allocate_check(ierr,'Br polar')
            allocate (Bth(n1,n2,n3),stat=ierr)
            call allocate_check(ierr,'Bth polar')
            read(60,err=98,iomsg=errmsg)Bz
            read(60,err=98,iomsg=errmsg)Br
            read(60,err=98,iomsg=errmsg)Bth
            Bz=Bz/bcon
            Br=Br/bcon
            Bth=Bth/bcon
            dtype='P3'
        end select
      Case('CART')
        nz=n1
        nx=n2
        ny=n3
        z0=x1
        y0=x2
        x0=x3
        dz=d1
        dx=d2
        dy=d3
        write(6,'("Axial magnetic field from file ",a132)')actual_name
!        write(6,'("z",i5,2f8.2)')nz,z0,dz
!        write(6,'("x",i5,2f8.2)')nx,x0,dx
!        write(6,'("y",i5,2f8.2)')ny,y0,dy
      WRITE(6,4020),NX,NY,NZ,DX,DY,DZ,X0,Y0,Z0
4020  FORMAT(' Inflector Magnetic Field',' NX=',I4,' NY=',I4,' NZ=',I4,/,1X,' DX=',F7.4,' DY=',F7.4,' DZ=',F7.4,' XP0=',F7.3,' YP0=',F7.3,' ZP0=',F7.3)
        
        select case(ndim)
          case(1)
            allocate (Bz(n1,1,1),stat=ierr)
            call allocate_check(ierr,'Bz cart')
            n2=1
            n3=1
            read(60,err=98,iomsg=errmsg)Bz
            Bz=Bz/bcon
            dtype='C1'
          case(2)
            allocate (Bz(n1,n2,1),stat=ierr)
            call allocate_check(ierr,'Bz cart')
            allocate (Bx(n1,n2,1),stat=ierr)
            call allocate_check(ierr,'Bx cart')
            read(60,err=98,iomsg=errmsg)Bz
            read(60,err=98,iomsg=errmsg)Bx
            Bz=Bz/bcon
            Bx=Bx/bcon
            dtype='C2'          
           case(3)
            allocate (Bz(n1,n2,n3),stat=ierr)
            call allocate_check(ierr,'Bz cart')
            allocate (Bx(n1,n2,n3),stat=ierr)
            call allocate_check(ierr,'Bx cart')
            allocate (By(n1,n2,n3),stat=ierr)
            call allocate_check(ierr,'By cart')
            read(60,err=98,iomsg=errmsg)Bz
            read(60,err=98,iomsg=errmsg)Bx
            read(60,err=98,iomsg=errmsg)By
            !write(77,'(f10.2,3f17.5)')(z0+(i-1)*dz,bz(i,1,1),bx(i,1,1),by(i,1,1),i=1,nz)
            Bz=Bz/bcon
            Bx=Bx/bcon
            By=By/bcon
            !write(77,'(2f10.5)')(z0+(i-1)*dz,bz(i,41,41),i=1,nz)
            dtype='C3'
        end select    
      case default
        write(6,'("Illegal field specification ",a5)')ftype
        stop '*** axial field file issue ***'
    end select
    close(60)
!    write(6,'("x,y,z,th= ",3(L1,1X),f10.5)')symx,symy,symz,sym_ang
    write(6,*)trim(comment)

    bLoaded=.true.
    return
    
99  continue
    write(6,'("Error opening axial bfield ")')
    write(6,'(A<csize>)')filename
97  write(6,'("Error reading line 1 of BAXIAL")')  
98  write(6,'(A<csize>)')errmsg
    return_msg=errmsg
    ierror=12
    bLoaded=.false.
    return
  end subroutine InputBaxial

! This is the main routine for computing the axial B field

  subroutine getBaxial(Y,bbx,bby,bbz,*)
    real(8),intent(in) :: Y(7) ! this is the integration vector, contains T,x,y,z
    real(8),intent(out):: bbx,bby,bbz ! field components divided by bcon
    ! local variables    
    real :: cofr(8),cofz(8),coft(8),cofx(8),cofy(8)
    real(8) ::sum1,sum2,sum3,ff,z,r,th,xi,yi,x,yy
    real :: bbr,bbt,b1,b2
    real :: thp ! angle after adjusting for field symmetry
    real :: sgnz,sgnx,sgny
    integer :: iz,ir,ith,ix,iy ! locations in the grid
    real :: fx,fy,fz,fr,fth ! fraction of a grid step
    integer :: i,j,k ! loop counters
    
    !if (.Not. bLoaded) call inputBaxial
    if (.Not. bLoaded) return 1
    select case (dtype)
      case ('P1', 'C1')
        if(symz) then
          z=abs(y(4))
          sgnz=sign(1.0D+0,y(4))
        else
          z=y(4)
          sgnz=1.0D+0
        end if
        r=sqrt(y(2)**2+y(3)**2)
        fz=(z-z0)/dz
        iz=fz
        fz=fz-iz
        if(iz+3 .gt. nz ) go to 999
	    if(iz .lt. 0) go to 999
        call intcf(fz,cofz)
        sum1=zero
        sum2=zero
        do i=1,4
          sum1=sum1+cofz(i)*bz(iz+i,1,1)
          sum2=sum2+cofz(i+4)*bz(iz+i,1,1)
        end do
        bbz=sum1
        bbr=-r*sum2*sgnz/dz/2.0d+0
        if(r .ne. 0.0) then
          bbx=bbr*y(2)/r
          bby=bbr*y(3)/r
        else
          bbx=zero
          bby=zero
        end if      
      case ('P2')
        r=sqrt(y(2)**2+y(3)**2)
        if(symz) then
          z=abs(y(4))
          sgnz=sign(1.0D+0,y(4))
        else
          z=y(4)
          sgnz=1.0D+0
        end if
        fz=(z-z0)/dz
        iz=fz
        fz=fz-iz
        if(iz+3 .gt. nz )go to 999
	    if(iz .lt. 0) go to 999
        fr=(r-r0)/dr
        ir=fr
        fr=fr-ir
        if(ir+3 .gt. nr ) go to 999
	    if(ir .lt. 0) go to 999
        call intcf(fz,cofz)
        call intcf(fr,cofr)
        sum1=0.0d+0
        sum2=0.0d+0
        do j=1,4
        do i=1,4
          ff=cofz(i)*cofr(j)
          sum1=sum1+ff*bz(iz+i-1,ir+j,1)
          sum2=sum2+ff*br(iz+i-1,ir+j,1)
        end do
        end do
        bbz=sum1
        bbr=sum2*sign(1.0d+0,y(4))
        if(r .ne. 0.0) then
          bbx=bbr*y(2)/r
          bby=bbr*y(3)/r
        else
          bbx=0.
          bby=0.
        end if      
      case ('P3')
        r=sqrt(y(2)**2+y(3)**2)
        if(symz) then
          z=abs(y(4))
          sgnz=sign(1.0D+0,y(4))
        else
          z=y(4)
          sgnz=1.0D+0
        end if
        th=atan2(y(3),y(2))-thr_rot
        if (th .lt. 1) th = th+ tpi
        thp=mod(th-th0,sym_ang)
        fz=(z-z0)/dz
        iz=fz
        fz=fz-iz
        if(iz+3 .gt. nz ) go to 999
	    if(iz .lt. 1) go to 999
        fr=(r-r0)/dr
        ir=fr
        fr=fr-ir
        if(ir+3 .gt. nr ) go to 999
	    if(ir .lt. 1) go to 999
        fth=thp/dth
        ith=fth
        fth=fth-ith
        if(ith+3 .gt. nth ) go to 999
	    if(ith .lt. 1) go to 999
        call intcf(fz,cofz)
        call intcf(fr,cofr)
        call intcf(fth,coft)
        sum1=0.0d+0
        sum2=0.0d+0
        sum3=0.0d+0
        do k=1,4
        do j=1,4
        do i=1,4
          ff=cofz(i)*cofr(j)*coft(k)
          sum1=sum1+ff*bz(iz+i-1,ir+j-1,ith+k-1)
          sum2=sum2+ff*br(iz+i-1,ir+j-1,ith+k-1)
          sum3=sum3+ff*bth(iz+i-1,ir+j-1,ith+k-1)
        end do
        end do
        end do
        bbz=sum1
              write(77,'(6f10.5)')r,th,z,sum1,sum2,sum3
        bbr=sum2*sign(1.0d+0,y(4))
        bbt=sum3*sign(1.0d+0,y(4))
        bbx=bbr*cos(th)-bbt*sin(th)
        bby=bbr*sin(th)+bbt*cos(th)      
      case ('C2')
        write(6,'("case not yet supported")')
        ierror=11
        return_msg='Case not yet supported'
        return
!        stop
      case('C3')
      ! x coord
      if(symy) then
        x=abs(y(2))
        sgnx=sign(1.0D+0,y(2))
      else
        x=y(2)
        sgnx=1.0D+0
      end if
      ! y coord
      if(symx) then
        yy=abs(y(3))
        sgny=sign(1.0D+0,y(3))
      else
        yy=y(3)
        sgny=1.0D+0
      end if
      ! z coord
      if(symz) then
        z=abs(y(4))
        sgnz=sign(1.0D+0,y(4))
      else
        z=y(4)
        sgnz=1.0D+0
      end if
      ! find interp param
      fz=(z-z0)/dz
      iz=fz
      fz=fz-iz
      if(iz+3 .gt. nz ) go to 999
	  if(iz .lt. 1) go to 999
      xi = x ! x inflector
      yi = yy ! y inflector
	  if(rotate) then
	    x = xi*cos_rot-yi*sin_rot !x map
	    yy= xi*sin_rot+yi*cos_rot !y map
	  endif
      fx=(x-x0)/dx
      ix=fx
      fx=fx-ix
      if(ix+3 .gt. nx ) go to 999
	  if(ix .lt. 1) go to 999
      fy=(yy-y0)/dy
      iy=fy
      fy=fy-iy
      if(iy+3 .gt. ny ) go to 999
	  if(iy .lt. 1) go to 999
      call intcf(fz,cofz)
      call intcf(fy,cofy)
      call intcf(fx,cofx)
        !write(77,'(3i5,4f10.5)')ix,iy,iz,cofz(1:4)
      sum1=zero
      sum2=zero
      sum3=zero
      do k=1,4
      do j=1,4
      do i=1,4
        ff=cofz(i)*cofx(j)*cofy(k)
        sum1=sum1+ff*bz(iz+i-1,ix+j-1,iy+k-1)
        sum2=sum2+ff*bx(iz+i-1,ix+j-1,iy+k-1)
        sum3=sum3+ff*by(iz+i-1,ix+j-1,iy+k-1)
      end do
      end do
      end do
        !write(77,'("x",3f10.5,4e16.8)')x,yy,z,(bz(iz+i,ix+1,iy+1),i=0,3)
      bbz=sum1
      bbx=sum2*sgnz*sgnx
      bby=sum3*sgnz*sgny
      if(rotate) then
        b1=bbx
        b2=bby
        bbx=b1*cos_rot+b2*sin_rot
        bby=-b1*sin_rot+b2*cos_rot
      endif
      ! zero the on axis values of Bx and By to stop inflector entrance from rotating
      if (abs(xi) .lt. 1e-4 .and. abs(yi) .lt. 1e-4) then
        bbx=0.0
        bby=0.0
      endif
    case default
      write(6,'("Illegal field type ")')
      stop 'Bfield error'     
    end select
    if(io_control80.active) write(80,'(9F12.5)')tcon*y(1),y(2:4),b0*bbx,b0*bby,b0*bbz
    return
    
999   write(6,'('' **** Off Field (BFIELD) **** '')')
      write(6,'('' iz,z,nz = '',i5,e13.5,i5,f10.5,1x,l1)')iz,z,nz,y(4),symz
      select case (dtype)
        case('P2', 'P3')
          write(6,'('' ir,r,nr = '',i5,f10.5,i5)')ir,r,nr
          write(6,'('' it,th,nth,thp = '',i5,f10.5,i5,f10.5)')ith,th*tcon,nth,thp*tcon
        case('C2','C3')
          write(6,'('' ix,x,nx = '',i5,e13.5,i5,f10.5,1x,l1)')ix,x,nx,y(2),symx
          write(6,'('' iy,y,ny = '',i5,e13.5,i5,f10.5,1x,l1)')iy,yy,ny,y(3),symy
      end select
      write(6,'(" cos,sin,xi,yi = ",l1,4f10.5)')rotate,cos_rot,sin_rot,xi,yi
      ierror=10
      return_msg=' **** Baxial Error Stop **** '
      return 1
      !stop ' **** Error Stop **** '    
  end subroutine getBaxial
  
! Called by Parameter_load to setup axial field read
  logical function SetBaxial
    real,external :: rdata_in
  
    th_rot=rdata_in(.false.,-360.0,360.,0.0)
    if(th_rot .ne. 0.0) then
      cos_rot=cosd(th_rot)
      sin_rot=sind(th_rot)
      thr_rot=th_rot/tcon
      rotate=.true.
    else
      thr_rot=0.0
      rotate=.false.
    endif
    call getLine(filename,*10,*10)
    !read(5,'(a<csize>)')filename ! get field file name
    filename=adjustl(filename)
    bLoaded=.false.
    SetBaxial=.false.
    return
10  SetBaxial=.true.
    write(6,'("Failure in SetBaxial to retrieve file name")')
  end function SetBaxial

end module