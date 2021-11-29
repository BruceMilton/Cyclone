
      PROGRAM CYCLONE
	use DFLIB
      use cyclone_data

      Character*200 :: run_name	! this determines the folder name for output
      Character*200 :: msg ! this is a return message  
      Integer*4 :: ierr ! error flag 0=okay      
	character*132 :: file_name    
c
c dos stuff
c
      run_name(1:1)=' '
      indx=nargs()
      if (indx .lt. 1) goto 5
      call getarg(1,run_name)
c
c first input line
c
 5    continue
      if (run_name(1:1) .eq. ' ') then
        write(6,'('' Enter run name '',$)')
        read(*,'(a32)')run_name
      endif
      if (len_trim(run_name)<1) goto 1900 ! blank entry = quit
      file_name=run_name(1:len_trim(run_name))//'.dat' 
      open(5,status='old',readonly,file=file_name,iostat=ios)
	IF(IOS .NE. 0) THEN
		WRITE(6,'('' Could not open file>'',A32)')run_name
	  run_name=' '
		go to 5
	endif
      
      ReadFile=.true.
      numInputLines=-1 ! means no data in array

      call cyclone_main(run_name,msg,ierr)
      
      if(ierr .ne. 0) write(6,'(A200)')msg

1900  STOP 'NORMAL EXIT'
      END