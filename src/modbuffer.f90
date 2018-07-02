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
module buffering

  use globvar

  implicit none

  contains
    subroutine dobuffercluster

      use globvar
      use ncdfpars
      use cellroutines

      implicit none
      
      include 'cdi.inc'      
      
      integer :: nIDs,globID
      real(kind=8), allocatable :: cl(:,:)
      
      ! data arrays
      real(kind=8), allocatable :: dat(:)          ! array for reading float from nc
      real(kind=8), allocatable :: dat2d(:,:)      ! array for doing the clustering
      
      bglobnIDs=0
      globID=1
      
      write(*,*)"======================================="
      write(*,*)"======== START CELL DETECTION ========="
      write(*,*)"======================================="
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to input file..."
    
      ! Get initial Information about grid and timesteps of both files
      CALL datainfo(bffile)
    
      ! Open the dataset 1
      streamID1=streamOpenRead(bffile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
    
      ! Set the variable IDs 1
      varID1=ivar
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)
    
      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(bfclfile)
      write(*,*)"---------"
    
      !! open new nc file for results
      ! define grid
      gridID2=gridCreate(GRID_GENERIC, nx*ny)
      CALL gridDefXsize(gridID2,nx)
      CALL gridDefYsize(gridID2,ny)
      CALL gridDefXvals(gridID2,xvals)
      CALL gridDefYvals(gridID2,yvals)
      CALL gridDefXunits(gridID2,TRIM(xunit))
      CALL gridDefYunits(gridID2,TRIM(yunit))
      zaxisID2=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID2, level)
      ! define variables
      vlistID2=vlistCreate()
      varID2=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID2,varID2,"cellID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each cell")
      CALL vlistDefVarUnits(vlistID2,varID2,"-")
      CALL vlistDefVarMissval(vlistID2,varID2,inmissval)
      CALL vlistDefVarDatatype(vlistID2,varID2,DATATYPE_INT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite(trim(bfclfile),FILETYPE_NC)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)
      
      ! Get data
      write(*,*)"======================================="
      write(*,*)"=== Find IDs for buffered cells ..."        

      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1 .OR. verbose)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if
    
        ! Allocate arrays for data storage
        allocate(dat(nx*ny))
    
        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)

        ! set time step for output file
        status=streamDefTimestep(streamID2,tsID)

        ! Read time step from input
        call streamReadVarSlice(streamID1,varID1,levelID,dat,nmiss1)
        
        ! cycle if field contains only missing values; but write it to output
        if(nmiss1==nx*ny)then
          nmiss2=nmiss1
          dat=inmissval
          CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
          deallocate(dat)
          cycle
        end if
        
        ! reshape array
        allocate(dat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        deallocate(dat)
    
        ! cluster the frame
        ! it is very important that at the end the cell IDs range from 1 to globnIDs without gaps
        allocate(cl(nx,ny))
        CALL clustering(dat2d,globID,globID,nIDs,cl,inmissval)
        !write(*,*)"clustering: Found ",nIDs," Cells"
        !write(*,*)nIDs,globnIDs,globID
        ! periodic boundaries
        if(nIDs>0 .AND. periodic)then
          globID=globID+1-nIDs
          CALL mergeboundarycells(cl,globID,globID,nIDs,inmissval)
          !write(*,*)"perbound: Found ",nIDs," Cells"
          !write(*,*)nIDs,globnIDs,globID
        end if
        if(nIDs.ne.0)globID=globID+1
        bglobnIDs=bglobnIDs+nIDs
        deallocate(dat2d)
        !write(*,*)nIDs,globnIDs,globID
        
        ! cells: reshape array and write to nc
        allocate(dat(nx*ny))
        CALL reshapeF2d(cl,nx,ny,dat)
        deallocate(cl)
        nmiss2=nmiss1
        CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
        deallocate(dat)

      end do
    
      ! close input and output
      CALL gridDestroy(gridID2)
      CALL vlistDestroy(vlistID2)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)
    
      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,*)"---------"
      write(*,*)"Cells found:",globnIDs
      write(*,*)"---------"
      
      if(globnIDs.eq.0)then
        write(*,*)"No cells found! We stop here!"
        write(*,*)"======================================="
        write(*,*)"======= FINISHED CELL DETECTION ======="
        write(*,*)"======================================="
        stop
      end if
    
      write(*,*)"======================================="
      write(*,*)"======= FINISHED CELL DETECTION ======="
      write(*,*)"======================================="
    
    end subroutine dobuffercluster

end module buffering
