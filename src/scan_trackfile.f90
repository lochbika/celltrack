!-------------------------------------------------------------------------------------------
!  ######  ######## ##       ##       ######## ########     ###     ######  ##    ## 
! ##    ## ##       ##       ##          ##    ##     ##   ## ##   ##    ## ##   ##  
! ##       ##       ##       ##          ##    ##     ##  ##   ##  ##       ##  ##   
! ##       ######   ##       ##          ##    ########  ##     ## ##       #####    
! ##       ##       ##       ##          ##    ##   ##   ######### ##       ##  ##   
! ##    ## ##       ##       ##          ##    ##    ##  ##     ## ##    ## ##   ##  
!  ######  ######## ######## ########    ##    ##     ## ##     ##  ######  ##    ##
!-------------------------------------------------------------------------------------------
! This file is part of celltrack
! Copyright: Kai Lochbihler (kai.lochbihler@knmi.nl)
!
! Any help will be appreciated :)
!
subroutine scan_trackfile(filename,ntracks,maxlen)

  use globvar, only : stdclen
  
  implicit none
  character(len=*), intent(in) :: filename
  character(len=stdclen) :: ttrack
  integer, intent(out) :: ntracks,maxlen
  integer :: length,riostat

  ! open the tracks file
  open(unit=1,file=filename,action="read",status="old")

  ! initialize counters
  maxlen=0
  length=0
  ntracks=0

  !now loop over all lines
  do
    ! read the line ;)
    read(1,'(A)',IOSTAT=riostat)ttrack
    ttrack=trim(ttrack)
    if(riostat.ne.0)exit
    !check if a new track begins
    if(ttrack(1:3)=="###")then
      read(1,*)
      ! buffer counter
      ntracks=ntracks+1
      if(maxlen<length)then
        maxlen=length
      end if
      length=0
    else
      length=length+1
    end if
  end do
  
  ! finally check one more time (for the last track) 
  ! whether it is longer than all other tracks..
  ! Credits go to Hylke de Vries (KNMI) for this fix
  if(maxlen<length)then
    maxlen=length
  end if

  if(maxlen<length)then
    maxlen=length
  end if


  close(unit=1)

end subroutine scan_trackfile
