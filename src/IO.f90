
      module IO
      implicit none
      private

! public routines
      PUBLIC :: input_io_control,input_io_control_simple,file_spacing,setRunDir,openFile,ioCleanup

! public variables

      TYPE,PUBLIC :: IO_Control
        logical :: Active =.FALSE.
        logical :: Formatted
        Integer :: Spacing ! gap between runs in file, -1 means EOF
        Integer :: Unit ! IO Unit number associated with this file
        Integer :: Lines ! Number of lines written in the current orbit
        Integer :: F34_State ! Indicator of file status that is written to unit 34
      END TYPE IO_Control
      
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control16
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control21
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control22
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control23
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control35
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control36
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control47
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control50
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control51
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control53
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control62
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control80
      TYPE(IO_Control),SAVE,PUBLIC :: IO_Control81

! Private variables
      character*200, private,save :: pathFile ! this will be the path plus the run_name
      character*200, private,save :: FileNames(99) !
      contains
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! called by parameter_load to setup the io control information
      
      Subroutine input_io_control(Control,iunit)

! dummy variables
        TYPE(IO_CONTROL) :: Control
        Integer, intent(in) :: iunit
! local variables
      integer,external :: idata_in,cdata_in
      logical,external :: ldata_in

      
          control%active=ldata_in(.false.,.true.)
	      control%formatted=ldata_in(.false.,.false.)
	      control%spacing=idata_in(.false.,-1,100,1)
	      control%unit=iunit

	      if(control.active) then
            if(control.formatted)then
  	          OPEN(control.unit,FORM='FORMATTED',STATUS='unknown',file=fileNames(control.unit))
  	        else
              OPEN(control.unit,FORM='UNFORMATTED',STATUS='unknown',file=fileNames(control.unit))	  
	        endif
	      else
	        close(control.unit)
	      endif
	      control.lines=0
	      if(control.active) then
	        if(control.formatted) then
	            control.F34_state=1
	        else
	            control.F34_state=-1
	        endif
	      else
	        control.f34_state = 0.0
	      endif
	      return
      end subroutine
      
     Subroutine input_io_control_simple(Control,iunit,state)

! dummy variables
        TYPE(IO_CONTROL) :: Control
        Integer, intent(in) :: iunit
        logical, intent(in) :: state
! local variables
      integer,external :: idata_in,cdata_in
      logical,external :: ldata_in

      
          control%active=state
	      control%formatted=.true.
	      control%spacing=0
	      control%unit=iunit

	      if(control.active) then
  	          OPEN(control.unit,FORM='FORMATTED',STATUS='unknown',file=fileNames(control.unit))
	      else
	        close(control.unit)
	      endif
	      control.lines=0
	      control.F34_state=1
	      return
      end subroutine      
      
      subroutine file_spacing
        
        call add_spaces(io_control16)
        call add_spaces(io_control35)
        call add_spaces(io_control36)
        call add_spaces(io_control47)                        
        call add_spaces(io_control50)
        call add_spaces(io_control51)
        call add_spaces(io_control53)
        call add_spaces(io_control62)
        call add_spaces(io_control80)

      end subroutine

! Add the defined inter run marker (EOF or empty lines)      
      subroutine add_spaces(control)
        ! Dummy variables
          TYPE(IO_CONTROL) :: Control
        ! local variables
          Integer :: i ! loop variable

        if (control.active) then
            if (control.formatted)  then
                if(control.spacing .eq. -1) then
                    endfile(control.unit)
                else
                    do i = 1, control.spacing
                        write(control.unit,'(" ")') 
                    enddo
                endif
            else
                endfile(control.unit)
            endif
            control.lines=0
        endif
      end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function setRunDir(run_name)
      
        use IFPORT
        
        character*(*) run_name
        integer :: iresult,i,iunit,irun,len,status
        character*200 path,msg,fileLoc
        character*20 ext
        logical :: exists
        
        path=FILE$CURDRIVE
        iresult=GETDRIVEDIRQQ(path)
        irun=len_trim(run_name)
        if(iresult .ne. 0) then
          pathFile=path(1:iresult)//'\'//run_name(1:irun)//'\'//run_name
        else
          write(6,'("Could not get current path")')
          pathFile=run_name
        endif

        do i =1,9
          write(fileNames(i),'("FORT",i1,".dat")',err=4)i
          fileNames(i)=run_name(1:irun)//'\'//fileNames(i)
4       enddo
        do i =10,99
          write(fileNames(i),'("FORT",i2,".dat")',err=5)i
          fileNames(i)=run_name(1:irun)//'\'//fileNames(i)
5       enddo
      ! look for a file called extensions.dat in the current folder
      ! if not there then look for the environment variable cyclone$dir

      fileLoc='extensions.dat'
      INQUIRE (FILE = fileLoc, EXIST = exists)
      IF (.NOT. exists) THEN
        call get_environment_variable ('cyclone$dir', fileLoc, len, status, .true.)
        if (status .ge. 2) then
          write (6,*) 'get_environment_variable failed: status = ', status
        end if
        if (status .eq. 1) then
          write (6,*) 'env var does not exist'
        end if
        if (status .eq. -1) then
          write (6,*) 'env var length = ', len, ' truncated to 200'
          len = 200
        end if
        if (status .eq. 0) then
          fileLoc=fileLoc(1:len)//'\'//'extensions.dat'
          exists=.true.
        endif
        if (status .eq. 0 .and. len .eq. 0) then
          write (6,*) 'env var exists  but has no value'
          exists=.false.
        end if
      endif      
!   if it exists then use the data in the file to assign files to units
      IF (exists) THEN
        write(6,'(A230)') 'Extensions file to be used is '// fileLoc
        open(95,file=fileLoc,readonly,err=11,iomsg=msg)
10      read(95,'(I5,A20)',err=11 ,end=12 ,iomsg=msg)iunit,ext
          if (iunit .le. 99 .and. iunit .ge. 1) then
            fileNames(iunit)=PathFile(1:len_trim(PathFile))//ext(1:len_trim(ext))
          else
            write(6,'("Bad unit number in ext.dat ",i5)')iunit
          endif
        go to 10
11      setRunDir=1
        write(6,*)'Extensions file not found'
        return
12      close(95)
      else
        ! no file found, so if FORT.6 is not defined then close so above messages are saved
        call get_environment_variable ('FORT.6', fileLoc, len, status, .true.)
        if(status .ne. 0) then
          close(6)
        endif   
      endif
      setRunDir=0
      end function
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function openFile(iunit,msg)
      
        integer,intent(in) :: iunit ! unit number of file to open
        character*(*),intent(out) :: msg
        character*200 test
        
        OPEN(iunit,FORM='FORMATTED',STATUS='unknown',file=fileNames(iunit),err=99,iomsg=msg)
        
        openFile=0
        return
99      openFile=1
        return
      end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! close all open files and any other housing keeping for files      
      integer function ioCleanup(dummy)
        integer,intent(in) :: dummy
        
        integer :: i
        do i=1,99
          close(i)
        enddo      
        iocleanup=dummy
      end function
      
      end module