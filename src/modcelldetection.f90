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
module celldetection
  
  use globvar

  implicit none
  
  
  contains 
    subroutine docelldetection
      
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
      
      globnIDs=0
      globID=1
      
      write(*,*)"======================================="
      write(*,*)"======== START CELL DETECTION ========="
      write(*,*)"======================================="
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to input file..."
    
      ! Get initial Information about grid and timesteps of both files
      CALL datainfo(ifile)
    
      ! Open the dataset 1
      streamID1=streamOpenRead(ifile)
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
    
      ! Get information about longitudes and latitudes and levels
      allocate(xvals(0:nx-1))
      allocate(yvals(0:ny-1))
      if(verbose)write(*,*)"Arrays for grid values successfully allocated!"
      nlev=zaxisInqSize(zaxisID1)
      !allocate(levels(1:nlev))
      !call zaxisInqLevels(zaxisID1,levels)
      nblon=gridInqXvals(gridID1,xvals)
      nblat=gridInqYvals(gridID1,yvals)
      diflon=(xvals(nx-1)-xvals(0))/(nblon-1)
      diflat=(yvals(ny-1)-yvals(0))/(nblat-1)
      level=zaxisInqLevel(zaxisID1,levelID)
      CALL vlistInqVarUnits(vlistID1,varID1,vunit)
      CALL gridInqXunits(gridID1,xunit)
      CALL gridInqYunits(gridID1,yunit)
      inmissval=vlistInqVarMissval(vlistID1,varID1)
      call vlistInqVarName(vlistID1,varID1,vname)
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
      
      write(*,*)"======================================="
      write(*,*)"=== INPUT SUMMARY:"
      write(*,*)"---------"
      write(*,*)"Input   : ",trim(ifile)
      write(*,*)"---------"
      write(*,'(A,1a12)')" VAR        : ",trim(vname)
      write(*,'(A,1a12)')" Unit       : ",trim(vunit)
      write(*,'(A,1f12.2)')" MissVal    : ",inmissval
      write(*,'(A,1i12)')" NX         : ",nx
      write(*,'(A,1f12.2)')" MIN X      : ",xvals(0)
      write(*,'(A,1f12.2)')" MAX X      : ",xvals(nblon-1)
      write(*,'(A,1f12.2)')" DIF X      : ",diflon
      write(*,'(A,1a12)')" Unit       : ",trim(xunit)
      write(*,'(A,1i12)')" NY         : ",ny
      write(*,'(A,1f12.2)')" MIN Y      : ",yvals(0)
      write(*,'(A,1f12.2)')" MAX Y      : ",yvals(nblat-1)
      write(*,'(A,1f12.2)')" DIF Y      : ",diflat
      write(*,'(A,1a12)')" Unit       : ",trim(yunit)
      write(*,'(A,1i8,1i8.6)')" START DATE : ",vdate(1),vtime(1)
      write(*,'(A,1i8,1I8.6)')" END DATE   : ",vdate(ntp),vtime(ntp)
      write(*,'(A,1i12)')" NTSTEPS    : ",ntp
      write(*,'(A,1i12)')" TSTEP      : ",tstep
      write(*,'(A,1i12)')" NLEV       : ",nlev
      write(*,'(A,1f12.2)')" SELLEV     : ",level
      write(*,*)"---------"
    
      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(outfile)
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
      streamID2=streamOpenWrite(trim(outfile),FILETYPE_NC)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)
      
      if(buffer>0)then

        write(*,*)"======================================="
        write(*,*)"=== CREATING OUTPUT ..."
        write(*,*)"Output  :     ",trim(bffile)
        write(*,*)"---------"
      
        !! open new nc file for results
        ! define grid
        gridID3=gridCreate(GRID_GENERIC, nx*ny)
        CALL gridDefXsize(gridID3,nx)
        CALL gridDefYsize(gridID3,ny)
        CALL gridDefXvals(gridID3,xvals)
        CALL gridDefYvals(gridID3,yvals)
        CALL gridDefXunits(gridID3,TRIM(xunit))
        CALL gridDefYunits(gridID3,TRIM(yunit))
        zaxisID3=zaxisCreate(ZAXIS_GENERIC, 1)
        CALL zaxisDefLevels(zaxisID3, level)
        ! define variables
        vlistID3=vlistCreate()
        varID3=vlistDefVar(vlistID3,gridID3,zaxisID3,TIME_VARIABLE)
        CALL vlistDefVarName(vlistID3,varID3,"bfarea")
        CALL vlistDefVarLongname(vlistID3,varID3,"buffered area around cells; 1=true, 0=false")
        CALL vlistDefVarUnits(vlistID3,varID3,"-")
        CALL vlistDefVarMissval(vlistID3,varID3,inmissval)
        CALL vlistDefVarDatatype(vlistID3,varID3,DATATYPE_INT32)
        ! copy time axis from input
        taxisID3=vlistInqTaxis(vlistID1)
        call vlistDefTaxis(vlistID3,taxisID3)
        ! Open the dataset for writing
        streamID3=streamOpenWrite(trim(bffile),FILETYPE_NC)
        if(streamID3<0)then
           write(*,*)cdiStringError(streamID3)
           stop
        end if
        ! Assign variables to dataset
        call streamDefVList(streamID3,vlistID3)
        
      end if
      
      ! save which time steps are completely filled with NA
      allocate(tsALLna(ntp))
      tsALLna=.false.

      ! Get data
      if(buffer==0)then
        write(*,*)"======================================="
        write(*,*)"=== Find continous cells ..."
      else
        write(*,*)"======================================="
        write(*,*)"=== Find continous cells and calculate buffers ..."        
      end if
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
        if(buffer>0)then
          status=streamDefTimestep(streamID3,tsID)
        end if

        ! Read time step from input
        call streamReadVarSlice(streamID1,varID1,levelID,dat,nmiss1)
        
        ! cycle if field contains only missing values; but write it to output
        if(nmiss1==nx*ny)then
          nmiss2=nmiss1
          dat=inmissval
          tsALLna(tsID+1)=.true.
          CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
          if(buffer>0)then
            nmiss3=nmiss1
            CALL streamWriteVar(streamID3,varID3,dat,nmiss3)
          end if
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
        ! delete small clusters/cells
        if(nIDs>0 .AND. minarea>0)then
          globID=globID+1-nIDs
          CALL delsmallcells(cl,globID,globID,nIDs,inmissval)
          !write(*,*)"dellsmall: Found ",nIDs," Cells"
          !write(*,*)nIDs,globnIDs,globID
        end if
        if(nIDs.ne.0)globID=globID+1
        globnIDs=globnIDs+nIDs
        deallocate(dat2d)
        !write(*,*)nIDs,globnIDs,globID
        
        ! buffer around remaining cells
        if(buffer>0)then
          allocate(bfmask(nx,ny))
          CALL distCellEdge(cl,bfmask)
        end if
        
        ! cells: reshape array and write to nc
        allocate(dat(nx*ny))
        CALL reshapeF2d(cl,nx,ny,dat)
        deallocate(cl)
        nmiss2=nmiss1
        CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
        deallocate(dat)

        ! cells: reshape array and write to nc
        if(buffer>0)then
          allocate(dat(nx*ny))
          CALL reshapeF2d(bfmask,nx,ny,dat)
          deallocate(bfmask)
          nmiss3=(nx*ny)-sum(dat)
          CALL streamWriteVar(streamID3,varID3,dat,nmiss3)
          deallocate(dat)
        end if
        
        
      end do
    
      ! close input and output
      CALL gridDestroy(gridID2)
      CALL vlistDestroy(vlistID2)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)
      if(buffer>0)then
        CALL gridDestroy(gridID3)
        CALL vlistDestroy(vlistID3)
        CALL streamClose(streamID3)
      end if

    
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
    
    end subroutine docelldetection
        
end module celldetection
