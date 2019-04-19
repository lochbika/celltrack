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

  use globvar, only : stdclen,stdfloattype

  implicit none

  ! Variables for netCDF actions
  integer :: gridID1,taxisID1,vlistID1,varID1,streamID1,tsID,zaxisID1
  integer :: gridID2,taxisID2,vlistID2,varID2,streamID2,zaxisID2,vuID,vvID
  integer :: gridID3,taxisID3,vlistID3,varID3,streamID3,zaxisID3,ssizeID
  integer :: gridID4,taxisID4,vlistID4,varID4,streamID4,zaxisID4
  integer :: nmiss1,nmiss2,nmiss3,nmiss4
  real(kind=stdfloattype) :: outmissval=-999.D0

  contains

    subroutine getVarInfo(infile)

      use globvar, only : ntp,status,tp,ivar,stdclen,vunit,vname,inmissval

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

      ! set name
      vname=ivar
      
      ! unit
      CALL vlistInqVarUnits(vlistID1,varID1,vunit)
      
      ! missing value
      inmissval=vlistInqVarMissval(vlistID1,varID1)
      if(NINT(inmissval).eq.0)inmissval=-123456789.D0
      
      call streamClose(streamID1)

    end subroutine getVarInfo
  
    subroutine getTaxisInfo(infile)

      use globvar, only : ntp,status,tp,ivar,stdclen,vdate,vtime

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
      taxisID1=vlistInqTaxis(vlistID1)

      ! get number of timesteps
      tp=0
      do tsID=0,999999
         status=streamInqTimestep(streamID1,tsID)
         if(status==0)exit
         tp=tp+1
      end do
      ntp=tp
      
      ! extract dates and times for all time steps
      allocate(vdate(ntp))
      allocate(vtime(ntp))
      
      do tsID=0,(ntp-1)
        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)
        ! read date and time
        vdate(tsID+1) = taxisInqVdate(taxisID1)
        vtime(tsID+1) = taxisInqVtime(taxisID1)
      end do
      
      call streamClose(streamID1)

    end subroutine getTaxisInfo
    
    subroutine gethorGrid(infile)

      use globvar, only : nx,ny,ivar,stdclen,ingrid,xvals,yvals,minx,maxx,miny,maxy,nblon,nblat,xunit,yunit

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

      ! Get information about number of grid points in each direction
      nx=gridInqXsize(gridID1)
      ny=gridInqYsize(gridID1)
      
      ! get type of grid
      ingrid=gridInqType(gridID1)

      ! extract grid point values (coordinates)
      if(ingrid.eq.GRID_GENERIC)then
        allocate(xvals(0:(nx-1)))
        allocate(yvals(0:(ny-1)))
      else if(ingrid.eq.GRID_CURVILINEAR)then
        allocate(xvals(0:(nx*ny-1)))
        allocate(yvals(0:(nx*ny-1)))
      else
        write(*,*)"ERROR: grid type not supported!"
        stop
      end if
      
      nblon=gridInqXvals(gridID1,xvals)
      nblat=gridInqYvals(gridID1,yvals)

      ! get minimum and maximum
      minx=MINVAL(xvals)
      maxx=MAXVAL(xvals)
      miny=MINVAL(yvals)
      maxy=MAXVAL(yvals)
      
      ! get units
      CALL gridInqXunits(gridID1,xunit)
      CALL gridInqYunits(gridID1,yunit)
      
      call streamClose(streamID1)

    end subroutine gethorGrid
    
    subroutine getZaxisInfo(infile)

      use globvar, only : ivar,stdclen,levels,level,levelID,nlev

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
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      nlev=zaxisInqSize(zaxisID1)
      allocate(levels(0:(nlev-1)))
      nlev=zaxisInqLevels(zaxisID1,levels)
      level=zaxisInqLevel(zaxisID1,levelID)
      
      call streamClose(streamID1)

    end subroutine getZaxisInfo
    
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
