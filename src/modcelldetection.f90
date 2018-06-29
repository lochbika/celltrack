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
    
    subroutine clustering(data2d,startID,finID,numIDs,tcl,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(out) :: finID,numIDs
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2)
      real(kind=8), intent(in) :: data2d(nx,ny),missval
      real(kind=8),intent(out) :: tcl(nx,ny)
      logical :: mask(nx,ny)
    
      ! initialize variables and arrays
      tcl=missval
      mask=.false.
      clID=startID
      numIDs=0
    
      ! mask values higher than threshold and if not missing value
      do y=1,ny
        do x=1,nx
          if(data2d(x,y)>thres .AND. data2d(x,y).ne.missval)then
          mask(x,y)=.true.
          end if
        end do
      end do
      
      ! check if there are any gridpoints to cluster
      if(ANY(mask))then
      
        ! assign IDs to continous cells
        do y=1,ny
          do x=1,nx
            neighb=missval
            if(mask(x,y))then
              ! gather neighbouring IDs; left,up
              if(x.ne.1)  neighb(1)=tcl((x-1),y)
              if(y.ne.1)  neighb(2)=tcl(x,(y-1))
              ! check if there is NO cluster around the current pixel; create new one
              if(ALL(neighb==missval))then
                tcl(x,y)=clID
                clID=clID+1
                numIDs=numIDs+1
              else
                ! both neighbours are in the same cluster
                if(neighb(1)==neighb(2).and.neighb(1).ne.missval)then
                  tcl(x,y)=neighb(1)
                end if
                ! both neighbors are in different clusters but none of them is (-999)
                if(neighb(1).ne.neighb(2) .and. neighb(1).ne.missval .and. neighb(2).ne.missval)then
                  numIDs=numIDs-1
                  tcl(x,y)=MINVAL(neighb)
                  ! update the existing higher cluster with the lowest neighbour
                  do cony=1,ny
                    do conx=1,nx
                      if(tcl(conx,cony)==MAXVAL(neighb))tcl(conx,cony)=MINVAL(neighb)
                    end do
                  end do
                end if
                ! both neighbors are in different clusters but ONE of them is empty(-999)
                if(neighb(1).ne.neighb(2) .and. (neighb(1)==missval .or. neighb(2)==missval))then
                  tcl(x,y)=MAXVAL(neighb)
                end if
              end if
            end if
          end do
        end do
        
        ! gather IDs and rename to gapless ascending IDs
        if(numIDs>0)then
          allocate(allIDs(numIDs))
          allIDs=missval
          clID=startID-1
          tp=1
          do y=1,ny
            do x=1,nx
              if(.NOT.ANY(allIDs==tcl(x,y)) .AND. tcl(x,y).ne.missval)then
                allIDs(tp)=tcl(x,y)
                tp=tp+1
              end if
            end do
          end do
      
          do i=1,tp-1
            clID=clID+1
            do y=1,ny
              do x=1,nx
                if(tcl(x,y)==allIDs(i))then
                  tcl(x,y)=clID
                end if
              end do
            end do
          end do
          deallocate(allIDs)
        end if
        
      end if
      ! return final cluster ID
      finID=clID
    end subroutine clustering

    subroutine mergeboundarycells(data2d,startID,finID,numIDs,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(inout) :: numIDs
      integer, intent(out) :: finID
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2)
      real(kind=8), intent(in) :: missval
      real(kind=8),intent(inout) :: data2d(nx,ny)
    
      ! initialize variables and arrays
    
      ! check in y direction if there are cells connected beyond boundaries
      do y=1,ny
        neighb=missval
        ! gather neighbouring IDs
        if(data2d(nx,y).ne.missval .AND. data2d(1,y).ne.missval)then
          neighb(1)=data2d(nx,y) ! this is the ID of the current gridpoint
          neighb(2)=data2d(1,y)  ! this is the ID of the neighbouring point over the edge of x
          if(neighb(1).ne.neighb(2))then
            ! now iterate the complete matrix and rename higher cluster ID to the lower ID
            do cony=1,ny
              do conx=1,nx
                if(data2d(conx,cony)==MAXVAL(neighb))data2d(conx,cony)=MINVAL(neighb)
              end do
            end do
            if(verbose)write(*,*)"y Cluster # ",MAXVAL(neighb)," replaced  by",MINVAL(neighb)
            ! one cluster was deleted
            numIDs=numIDs-1            
          end if
        end if
      end do
      ! check in x direction if there are cells connected beyond boundaries
      do x=1,nx
        neighb=missval
        ! gather neighbouring IDs
        if(data2d(x,ny).ne.missval .AND. data2d(x,1).ne.missval)then
          neighb(1)=data2d(x,ny) ! this is the ID of the current gridpoint
          neighb(2)=data2d(x,1)  ! this is the ID of the neighbouring point over the edge of y
          if(neighb(1).ne.neighb(2))then
            ! now iterate the complete matrix and rename higher cluster ID to the lower ID
            do cony=1,ny
              do conx=1,nx
                if(data2d(conx,cony)==MAXVAL(neighb))data2d(conx,cony)=MINVAL(neighb)
              end do
            end do
            if(verbose)write(*,*)"x Cluster # ",MAXVAL(neighb)," replaced  by",MINVAL(neighb)
            ! one cluster was deleted
            numIDs=numIDs-1            
          end if
        end if
      end do
      
      ! gather IDs and rename to gapless ascending IDs
      allocate(allIDs(numIDs))
      allIDs=-999
      tp=1
      do y=1,ny
        do x=1,nx
          if(.NOT.ANY(allIDs==data2d(x,y)) .AND. data2d(x,y).ne.missval)then
            allIDs(tp)=data2d(x,y)
            tp=tp+1
          end if
        end do
      end do
  
      clID=startID-1
      do i=1,tp-1
        clID=clID+1
        do y=1,ny
          do x=1,nx
            if(data2d(x,y)==allIDs(i))then
              data2d(x,y)=clID
            end if
          end do
        end do
      end do
      deallocate(allIDs)
        
      ! return final cluster ID
      finID=clID
    end subroutine mergeboundarycells
    
    subroutine delsmallcells(data2d,startID,finID,numIDs,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres,minarea
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(inout) :: numIDs
      integer, intent(out) :: finID
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2),area(numIDs)
      real(kind=8), intent(in) :: missval
      real(kind=8),intent(inout) :: data2d(nx,ny)
    
      ! calculate each cells area
      area=0
      do x=1,nx
        do y=1,ny
          if(data2d(x,y).ne.missval)then
          !write(*,*)data2d(x,y),x,y
          area(INT(data2d(x,y))+1-startID) = area(INT(data2d(x,y))+1-startID) + 1
          end if
        end do
      end do
      
      ! how many clusters will be deleted?
      do i=1,numIDs
        if(area(i)<minarea)numIDs=numIDs-1 
      end do
      
      ! delete clusters with area smaller than minarea
      do y=1,ny
        do x=1,nx
          if(data2d(x,y).ne.missval)then
            if(area(INT(data2d(x,y))+1-startID)<minarea)then
              data2d(x,y)=missval
              if(verbose)write(*,*)"Cluster # ",i," deleted"
            end if
          end if
        end do
      end do      
      
      if(numIDs>0)then ! otherwise all clusters were deleted!
        ! gather IDs and rename to gapless ascending IDs
        allocate(allIDs(numIDs))
        allIDs=missval
        tp=1
        do y=1,ny
          do x=1,nx
            if(.NOT.ANY(allIDs==data2d(x,y)) .AND. data2d(x,y).ne.missval)then
              allIDs(tp)=data2d(x,y)
              tp=tp+1
            end if
          end do
        end do
    
        clID=startID-1
        do i=1,tp-1
          clID=clID+1
          do y=1,ny
            do x=1,nx
              if(data2d(x,y)==allIDs(i))then
                data2d(x,y)=clID
              end if
            end do
          end do
        end do
        deallocate(allIDs)
        
        ! return final cluster ID
        finID=clID
      else
        finID=startID
      end if
    
    end subroutine delsmallcells
    
    subroutine distCellEdge(data2d,mask)

      use globvar, only : periodic,y,x,i,tp,nx,ny,buffer

      implicit none
      ! input
      real(kind=8), intent(in) :: data2d(nx,ny) ! the input 2D field; we expect the field with cellIDs
      ! output
      real(kind=8), intent(out):: mask(nx,ny)   ! output: integer field; 1 if grid point is in buffer, else 0
      ! logical masks
      logical                  :: prcmsk(nx,ny) ! logical field; true if value > thres
      logical                  :: edgmsk(nx,ny) ! logical field; true if grid point is cell edge
      ! etc.
      integer                  :: nedg          ! number of points with edges and prc
      integer,allocatable      :: coorde(:,:)   ! coordinates of edge points
      real(kind=8)             :: dist,mindist,distxo,distxp
      real(kind=8)             :: distyo,distyp,distx,disty
      real(kind=8)             :: tcl(nx,ny)    ! field holding the distances
          
      ! initialize variables and arrays
      tcl=0.0
      prcmsk=.false.
      edgmsk=.false.
      mask=0
    
      ! mask values higher than threshold
      do y=1,ny
        do x=1,nx
          if(data2d(x,y)>0)prcmsk(x,y)=.true.
        end do
      end do
    
      ! mask cell edges
      do y=1,ny
        do x=1,nx
          if(prcmsk(x,y))then
            if(x.ne.1)then
              if(.NOT.prcmsk(x-1,y))edgmsk(x,y)=.true.
            else
              edgmsk(x,y)=.true.
            end if
            if(y.ne.1)then
              if(.NOT.prcmsk(x,y-1))edgmsk(x,y)=.true.
            else
              edgmsk(x,y)=.true.
            end if
            if(x.ne.nx)then
              if(.NOT.prcmsk(x+1,y))edgmsk(x,y)=.true.
            else
              edgmsk(x,y)=.true.
            end if
            if(y.ne.ny)then
              if(.NOT.prcmsk(x,y+1))edgmsk(x,y)=.true.
            else
              edgmsk(x,y)=.true.
            end if
          end if
        end do
      end do
    
      ! accound for periodic boundaries
      if(periodic)then
        do x=1,nx
          if(prcmsk(x,1))then
            if(.NOT.prcmsk(x,ny))edgmsk(x,1)=.true.
          end if
          if(prcmsk(x,ny))then
            if(.NOT.prcmsk(x,1))edgmsk(x,ny)=.true.
          end if
        end do
        do y=1,ny
          if(prcmsk(1,y))then
            if(.NOT.prcmsk(nx,y))edgmsk(1,y)=.true.
          end if
          if(prcmsk(nx,y))then
            if(.NOT.prcmsk(1,y))edgmsk(nx,y)=.true.
          end if
        end do
      end if
    
      ! count number of edge points
      nedg=0
      do y=1,ny
        do x=1,nx
          if(edgmsk(x,y))nedg=nedg+1
        end do
      end do
    
      ! create vectorized version of the coordinates of edge points
      allocate(coorde(nedg,2))
      tp=1
      do y=1,ny
        do x=1,nx
          if(edgmsk(x,y))then
            coorde(tp,1)=x
            coorde(tp,2)=y
            tp=tp+1
          end if
        end do
      end do
    
      ! find the closest edge point for each grid point
      do y=1,ny
        do x=1,nx
          if(.NOT.edgmsk(x,y))then
            ! calculate distances to edge points and find shortest
            mindist=huge(mindist)
            do i=1,nedg
              if(periodic)then
                ! ok, now it's getting tricky because we don't know whether selCL
                ! crossed the boundaries to become clID.
                ! At this step I calculate two kinds of distances between the two cells
                ! 1. the direct distance
                !    this assumes selCL did not cross the boundaries
                ! 2. the distance between clID and the projection of selCL
                !    this assumes that selCL crossed the boundaries
    
                ! clID=grid point
                ! selCL=edge point
    
                distxo=x-coorde(i,1) ! 1 in x direction
                distyo=y-coorde(i,2) ! 1 in y direction
                ! 2 in x direction
                if(x>=coorde(i,1))then
                  distxp=x-(coorde(i,1)+nx)
                else
                  distxp=(x+nx)-coorde(i,1)
                end if
                ! 2 in y direction
                if(y>=coorde(i,2))then
                  distyp=y-(coorde(i,2)+ny)
                else
                  distyp=(y+ny)-coorde(i,2)
                end if
                if( abs(distxo) .gt. abs(distxp) )then
                  ! this means the cell crossed the boundaries
                  distx=distxp
                else
                  ! no boundary crossing: just do the normal calculation
                  distx=distxo
                end if
                if( abs(distyo) .gt. abs(distyp) )then
                  ! this means the cell crossed the boundaries
                  disty=distyp
                else
                  ! no boundary crossing: just do the normal calculation
                  disty=distyo
                end if
              else
                distx=x-coorde(i,1)
                disty=y-coorde(i,2)
              end if
              ! calculate the final distance between this pair
              dist = sqrt(distx**2 + disty**2)
              if(dist<mindist)mindist=dist
              ! if the distance is 1 we found our closest point; shorter than 1 is not possible
              if(mindist.eq.1)exit
            end do
            if(prcmsk(x,y))then
              tcl(x,y)=mindist*(-1)
            else
              tcl(x,y)=mindist
            end if
          else
            tcl(x,y)=0.0D0
          end if
        end do
      end do

      deallocate(coorde)
      
      ! mask points within a certain distance
      do y=1,ny
        do x=1,nx
          if(tcl(x,y)<=buffer)then
            mask(x,y)=1
          end if
        end do
      end do
        
    end subroutine distCellEdge

end module celldetection
