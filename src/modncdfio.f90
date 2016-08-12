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
  implicit none

  ! Variables for netCDF actions
  integer :: gridID1,taxisID1,vlistID1,varID1,streamID1,tsID,zaxisID1
  integer :: gridID2,taxisID2,vlistID2,varID2,streamID2,zaxisID2,vuID,vvID
  integer :: gridID3,taxisID3,vlistID3,varID3,streamID3,zaxisID3,ssizeID
  integer :: nmiss1,nmiss2,nmiss3
  real(kind=8) :: inmissval,outmissval
  character(len=1000) :: vunit,xunit,yunit,vname

  contains
    subroutine datainfo(ifile)

      use globvar, only : nx,ny,ntp,status,tp

      implicit none

      include 'cdi.inc'

      integer :: gridID,taxisID,varID,streamID,tsID,vlistID
      character(len=500), intent(in) :: ifile

      ! Open the dataset
      streamID=streamOpenRead(ifile)
      if(streamID<0)then
         write(0,*)cdiStringError(streamID)
         stop
      end if
      ! Set the variable IDs
      varID=0
      vlistID=streamInqVlist(streamID)
      gridID=vlistInqVarGrid(vlistID,varID)
      taxisID=vlistInqTaxis(vlistID)
      ! Get information about number of time steps, longitudes and latitudes
      nx=gridInqXsize(gridID)
      ny=gridInqYsize(gridID)
      ! get number of timesteps
      tp=0
      do tsID=0,999999
         status=streamInqTimestep(streamID,tsID)
         if(status==0)exit
         tp=tp+1
      end do
      ntp=tp
      call streamClose(streamID)

    end subroutine datainfo

end module ncdfpars
