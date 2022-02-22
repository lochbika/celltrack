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
      write(*,*)"=== START BUFFERED CELL DETECTION ====="
      write(*,*)"======================================="
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to input file..."
    
      ! Get initial Information about grid and timesteps of both files
      CALL getVarInfo(bffile)
    
      ! Open the dataset 1
      streamID1=streamOpenRead(bffile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
    
      ! Set the variable IDs 1
      varID1=getVarIDbyName(bffile,"bfarea")
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
      CALL vlistDefVarDatatype(vlistID2,varID2,CDI_DATATYPE_INT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite(trim(bfclfile),CDI_FILETYPE_NC)
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
        ! it is very important that at the end the cell IDs range from 1 to bglobnIDs without gaps
        allocate(cl(nx,ny))
        CALL clustering(dat2d,globID,nIDs,cl,inmissval)
        !write(*,*)"clustering: Found ",nIDs," Cells"
        !write(*,*)nIDs,globnIDs,globID
        ! periodic boundaries
        if(nIDs>0 .AND. periodic)then
          globID=globID+1-nIDs
          CALL mergeboundarycells(cl,globID,nIDs,outmissval)
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
    
      write(*,*)"======================================="
      write(*,*)"== FINISHED BUFFERED CELL DETECTION ==="
      write(*,*)"======================================="
    
    end subroutine dobuffercluster

    subroutine linkbfcells

      use ncdfpars
      use globvar

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=8), allocatable :: dat(:),pdat(:)          ! array for reading float from nc
      real(kind=8), allocatable :: dat2d(:,:),pdat2d(:,:)  ! array for doing the clustering

      real(kind=8), allocatable :: wsum(:)

      write(*,*)"======================================="
      write(*,*)"=========== CELL STATISTICS ==========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== collecting timesteps ..."
      write(*,*)"---------"

      !!!!!!!!!
      ! number of unique IDs is globnIDs
      ! loop time steps to find all IDs and their timesteps

      ! Ids range from 1 to globnIDs
      allocate(bclIDs(bglobnIDs))
      bclIDs=-1
      do i=1,bglobnIDs
        bclIDs(i)=i
      end do

      ! Open the buffered cells file
      streamID2=streamOpenRead(bfclfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=getVarIDbyName(outfile,"cellID")
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      allocate(btsclID(bglobnIDs))
      allocate(bcldate(bglobnIDs))
      allocate(bcltime(bglobnIDs))
      btsclID=-1
      tp=1
      ltp=1
      bcldate=-1
      bcltime=-1

      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if

        ! Allocate arrays for data storage
        allocate(dat(nx*ny))

        ! Set time step for input file
        status=streamInqTimestep(streamID2,tsID)

        ! Read time step from input
        call streamReadVar(streamID2,varID2,dat,nmiss2)

        !cycle if all values are missing
        if(nmiss2==nx*ny)then
          if(verbose)write(*,*)"NO clusters in timestep:  ",tsID+1
          deallocate(dat)
          cycle
        end if

        ! assign time steps
        do i=1,nx*ny
          if(dat(i).ne.inmissval)then
            if(btsclID(INT(dat(i)))==-1)then
              btsclID(INT(dat(i)))=tsID+1
              bcldate(INT(dat(i)))=vdate(tsID+1)
              bcltime(INT(dat(i)))=vtime(tsID+1)
            end if
            if(verbose)write(*,*)"cluster: ",bclIDs(INT(dat(i)))," is at timestep: ",btsclID(INT(dat(i)))
            if(verbose)write(*,*)bcldate(INT(dat(i))),bcltime(INT(dat(i)))
          end if
        end do
        deallocate(dat)
      end do

      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"=== calculating area, (weighted) center of mass, peak and average values ..."
      write(*,*)"---------"

      ! Open the buffered cells file
      streamID2=streamOpenRead(bfclfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=getVarIDbyName(bfclfile,"cellID")
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      ! Open the cells file
      streamID1=streamOpenRead(outfile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
      varID1=getVarIDbyName(outfile,"cellID")
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      ! allocate and initialize arrays for cell statistics
      allocate(bclarea(bglobnIDs))
      allocate(btouchb(bglobnIDs))
      allocate(bclcmass(bglobnIDs,2))
      allocate(bwclcmass(bglobnIDs,2))
      allocate(bclavint(bglobnIDs))
      allocate(bclpint(bglobnIDs))
      allocate(wsum(bglobnIDs))
      clarea=0
      touchb=.false.
      clavint=0
      clpint=0
      wsum=0
      clcmass=0
      wclcmass=0

      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if

        ! if this time step has no values: set touchb to true for all cells in previous ts
        ! do this also for cells in next time step
        ! because it could be that this interrupts tracks in the middle of their life time
        if(tsALLna(tsID+1))then
          do i=1,globnIDs
            if(tsclID(i)==tsID .OR. tsclID(i)==tsID+2)touchb(i)=.true.
          end do
        end if

        ! cycle if there are no cells in this time step
        if(.NOT.ANY(tsclID==tsID+1))cycle

        ! Set time step for input files
        status=streamInqTimestep(streamID1,tsID)
        status=streamInqTimestep(streamID2,tsID)

        ! Allocate arrays for data storage
        allocate(dat(nx*ny))
        allocate(pdat(nx*ny))

        ! Read time step from input
        call streamReadVar(streamID2,varID2,dat,nmiss2)
        call streamReadVarSlice(streamID1,varID1,levelID,pdat,nmiss1)

        ! reashape to 2d; better for center of mass calculation
        allocate(dat2d(nx,ny),pdat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        CALL reshapeT2d(pdat,nx,ny,pdat2d)
        deallocate(dat,pdat)

        ! now loop dat2d and calculate statistics
        do y=1,ny
          do x=1,nx
            if(dat2d(x,y)==inmissval)cycle
            ! this cell touches missing values?
            if(x.ne.1)then
              if(pdat2d(x-1,y)==inmissval)touchb(INT(dat2d(x,y))) = .true.
            end if
            if(y.ne.1)then
              if(pdat2d(x,y-1)==inmissval)touchb(INT(dat2d(x,y))) = .true.
            end if
            if(x.ne.nx)then
              if(pdat2d(x+1,y)==inmissval)touchb(INT(dat2d(x,y))) = .true.
            end if
            if(y.ne.ny)then
              if(pdat2d(x,y+1)==inmissval)touchb(INT(dat2d(x,y))) = .true.
            end if
            ! this cell touches time or space boundaries?
            if(periodic)then
              if(tsID==0 .OR. tsID==ntp-1)touchb(INT(dat2d(x,y))) = .true.
            else
              if(y==1 .OR. x==1 .OR. x==nx .OR. y==ny .OR. tsID==0 .OR. tsID==ntp-1)touchb(INT(dat2d(x,y))) = .true.
            end if 
            ! cell area
            clarea(INT(dat2d(x,y))) = clarea(INT(dat2d(x,y))) + 1
            ! average intensity
            clavint(INT(dat2d(x,y))) = clavint(INT(dat2d(x,y))) + pdat2d(x,y)
            ! peak intensity
            if(clpint(INT(dat2d(x,y)))<pdat2d(x,y))clpint(INT(dat2d(x,y))) = pdat2d(x,y)
            ! centers of masses
            clcmass(INT(dat2d(x,y)),1) = clcmass(INT(dat2d(x,y)),1) + x
            clcmass(INT(dat2d(x,y)),2) = clcmass(INT(dat2d(x,y)),2) + y
            wclcmass(INT(dat2d(x,y)),1) = wclcmass(INT(dat2d(x,y)),1) + x * pdat2d(x,y)
            wclcmass(INT(dat2d(x,y)),2) = wclcmass(INT(dat2d(x,y)),2) + y * pdat2d(x,y)
            wsum(INT(dat2d(x,y))) = wsum(INT(dat2d(x,y))) + pdat2d(x,y)
          end do
        end do
        deallocate(dat2d,pdat2d)
      end do

      ! divide by area
      do i=1,globnIDs
        ! average intensity
        clavint(i)=clavint(i)/clarea(i)
        ! center of mass
        clcmass(i,1)=clcmass(i,1)/clarea(i)
        clcmass(i,2)=clcmass(i,2)/clarea(i)
        ! weighted center of mass
        wclcmass(i,1)=wclcmass(i,1)/wsum(i)
        wclcmass(i,2)=wclcmass(i,2)/wsum(i)
      end do

      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,'(A,1i12)')" AVerage cell area(gridpoints):",SUM(clarea)/globnIDs
      write(*,'(A,1i12)')" Minimum cell area(gridpoints):",MINVAL(clarea)
      write(*,'(A,1i12)')" Maximum cell area(gridpoints):",MAXVAL(clarea)
      write(*,'(A,1f12.6)')" TOtal average value          :",SUM(clavint)/globnIDs
      write(*,'(A,1f12.6)')" Mininmum avergae value       :",MINVAL(clavint)
      write(*,'(A,1f12.6)')" Maxinmum avergae value       :",MAXVAL(clavint)
      write(*,'(A,1f12.6)')" AVerage peak value           :",SUM(clpint)/globnIDs
      write(*,'(A,1f12.6)')" Mininmum peak value          :",MINVAL(clpint)
      write(*,'(A,1f12.6)')" Maxinmum peak value          :",MAXVAL(clpint)
      write(*,*)"---------"

      write(*,*)"======================================="
      write(*,*)"=== writing stats to file cell_stats.txt ..."
      write(*,*)"---------"

      ! write the stats to file
      open(unit=1,file="cell_stats.txt",action="write",status="replace")
      write(1,*)"     clID    tsclID    clarea       clcmassX"//&
      &"       clcmassY      wclcmassX      wclcmassY            peakVal"//&
      &"              avVal  touchb      date(YYYYMMDD)        time(hhmmss)"
      do i=1,globnIDs
        write(1,'(3i10,4f15.6,2f19.12,1L8,2i20.6)')clIDs(i),tsclID(i),clarea(i),clcmass(i,:), &
          & wclcmass(i,:),clpint(i),clavint(i),touchb(i),cldate(i),cltime(i)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"========= FINISHED STATISTICS ========="
      write(*,*)"======================================="

    end subroutine linkbfcells

end module buffering
