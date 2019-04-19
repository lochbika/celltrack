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
module ncdfpars

  use globvar, only : stdclen

  implicit none

  ! Variables for netCDF actions
  integer :: gridID1,taxisID1,vlistID1,varID1,streamID1,tsID,zaxisID1
  integer :: gridID2,taxisID2,vlistID2,varID2,streamID2,zaxisID2,vuID,vvID
  integer :: gridID3,taxisID3,vlistID3,varID3,streamID3,zaxisID3,ssizeID
  integer :: gridID4,taxisID4,vlistID4,varID4,streamID4,zaxisID4
  integer :: nmiss1,nmiss2,nmiss3,nmiss4
  real(kind=8) :: inmissval
  real(kind=8) :: outmissval=-999.D0
  character(len=stdclen) :: vunit,xunit,yunit,vname

  contains
    subroutine datainfo(infile)

      use globvar, only : nx,ny,ntp,status,tp,ivar,stdclen

      implicit none

      include 'cdi.inc'

      character(len=*), intent(in) :: infile

      ! Open the dataset
      streamID1=streamOpenRead(infile)
      if(streamID1<0)then
         write(0,*)cdiStringError(streamID1)
         stop
      end if
      ! Set the variable IDs
      varID1=getVarIDbyName(infile,ivar)
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      ! Get information about number of time steps, longitudes and latitudes
      nx=gridInqXsize(gridID1)
      ny=gridInqYsize(gridID1)
      ! get number of timesteps
      tp=0
      do tsID=0,999999
         status=streamInqTimestep(streamID1,tsID)
         if(status==0)exit
         tp=tp+1
      end do
      ntp=tp
      call streamClose(streamID1)

    end subroutine datainfo
    
    integer function getVarIDbyName(infile,name)
    
      ! this function returns the variable ID for a certain variable (name)

      implicit none

      include 'cdi.inc'

      character(len=*) , intent(in) :: infile
      character(len=*), intent(in) :: name
      character(len=stdclen) :: testname
      
      integer :: varID,streamID,vlistID,nVars,i

      ! Open the dataset
      streamID=streamOpenRead(infile)
      if(streamID<0)then
         write(0,*)cdiStringError(streamID)
         stop
      end if

      ! Get the total number of variables
      vlistID=streamInqVlist(streamID)
      nVars=vlistNvars(vlistID)
      
      ! iterate Variables and check
      do i=0,(nVars-1)
        ! extract variable name
        CALL vlistInqVarName(vlistID,i,testname)
        if(TRIM(testname).eq.TRIM(name))then
          varID=i
          exit
        end if
        if(i==(nVars-1))varID=-100
      end do

      call streamClose(streamID)
      
      if(varID.eq.-100)then
        write(*,*)"Error: Variable ",TRIM(name)," not found in data set!"
        stop
      end if
      
      getVarIDbyName=VarID

    end function getVarIDbyName

!     subroutine createOutputFile(infile,outfile)
!     
!       this subroutine creates the netCDF file (outfile) which holds all output variables
!       initially, it will create an netCDF file which is empty
!       variables can be added with a seperate subroutine
!       
!       grid, time axis and zaxis will be copied from (infile)
!     
!       use globvar, only : nx,ny,ntp,status,tp,ivar
! 
!       implicit none
! 
!       include 'cdi.inc'
! 
!       character(len=800), intent(in) :: infile
!       character(len=800), intent(in) :: outfile
! 
!       Open the dataset
!       streamID1=streamOpenRead(infile)
!       if(streamID1<0)then
!          write(0,*)cdiStringError(streamID1)
!          stop
!       end if
!       Set the variable IDs
!       varID1=ivar
!       vlistID1=streamInqVlist(streamID1)
!       gridID1=vlistInqVarGrid(vlistID1,varID1)
!       taxisID1=vlistInqTaxis(vlistID1)
!       zaxisID1=vlistInqVarZaxis(vlistID1,varID1)
! 
!       Get information about number of time steps, longitudes and latitudes
!       nx=gridInqXsize(gridID1)
!       ny=gridInqYsize(gridID1)
!       get number of timesteps
!       tp=0
!       do tsID=0,999999
!          status=streamInqTimestep(streamID1,tsID)
!          if(status==0)exit
!          tp=tp+1
!       end do
!       ntp=tp
!       call streamClose(streamID1)
! 
!     end subroutine createOutputFile
    
end module ncdfpars
