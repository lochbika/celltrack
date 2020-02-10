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
module cellstatistics

  use globvar

  implicit none

  contains
    subroutine calccellstatistics

      use ncdfpars
      use globvar

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=stdfloattype), allocatable :: dat(:),pdat(:)          ! array for reading float from nc
      real(kind=stdfloattype), allocatable :: dat2d(:,:),pdat2d(:,:)  ! array for doing the clustering

      real(kind=stdfloattype), allocatable :: coords(:,:,:)     ! array for holding all cells coordinates
      real(kind=stdfloattype), allocatable :: coordsgrd(:,:,:)  ! array for holding all cells coordinates (gridpoints)
      real(kind=stdfloattype), allocatable :: ints(:,:)         ! array for holding all cells intenities
      integer, allocatable :: cellcnt(:)             ! array for holding a counter for each cell

      logical :: clbnd ! true if a cell touches domain boundaries
      logical :: projx,projy ! shall we project points across a boundary
      integer :: nrupbndx,nrupbndy ! how many grid points are closer to the upper boundaries

      integer :: minclID,maxclID,maxarea

      real(kind=stdfloattype), allocatable :: wsum(:)

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
      allocate(clIDs(globnIDs))
      clIDs=-1
      do i=1,globnIDs
        clIDs(i)=i
      end do

      ! Open the cells file
      streamID2=streamOpenRead(outfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=getVarIDbyName(outfile,"cellID")
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      allocate(tsclID(globnIDs))
      allocate(cldate(globnIDs))
      allocate(cltime(globnIDs))
      tsclID=-1
      tp=1
      ltp=1
      cldate=-1
      cltime=-1

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
          if(dat(i).ne.outmissval)then
            if(tsclID(INT(dat(i)))==-1)then
              tsclID(INT(dat(i)))=tsID+1
              cldate(INT(dat(i)))=vdate(tsID+1)
              cltime(INT(dat(i)))=vtime(tsID+1)
            end if
            if(verbose)write(*,*)"cluster: ",clIDs(INT(dat(i)))," is at timestep: ",tsclID(INT(dat(i)))
            if(verbose)write(*,*)cldate(INT(dat(i))),cltime(INT(dat(i)))
          end if
        end do
        deallocate(dat)
      end do

      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"=== calculating area, peak and average values ..."
      write(*,*)"---------"

      !!!!!!!!!
      ! find center of mass and area of clusters

      ! Open the cells file
      streamID2=streamOpenRead(outfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=getVarIDbyName(outfile,"cellID")
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      ! Open the original data file
      streamID1=streamOpenRead(ifile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
      varID1=getVarIDbyName(ifile,ivar)
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      ! allocate and initialize arrays for cell statistics
      allocate(clarea(globnIDs))
      allocate(clareagrd(globnIDs))
      allocate(touchb(globnIDs))
      allocate(clcmass(globnIDs,2))
      allocate(wclcmass(globnIDs,2))
      allocate(clcmassgrd(globnIDs,2))
      allocate(wclcmassgrd(globnIDs,2))
      allocate(clavint(globnIDs))
      allocate(clpint(globnIDs))
      allocate(wsum(globnIDs))
      clarea=0.D0
      clareagrd=0
      touchb=.false.
      clavint=0
      clpint=0
      wsum=0
      clcmass=0
      wclcmass=0
      clcmassgrd=0
      wclcmassgrd=0

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

        ! reashape to 2d
        allocate(dat2d(nx,ny),pdat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        CALL reshapeT2d(pdat,nx,ny,pdat2d)
        deallocate(dat,pdat)

        ! now loop dat2d and calculate statistics
        do y=1,ny
          do x=1,nx
            if(dat2d(x,y)==outmissval)cycle
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
            clareagrd(INT(dat2d(x,y))) = clareagrd(INT(dat2d(x,y))) + 1
            clarea(INT(dat2d(x,y))) = clarea(INT(dat2d(x,y))) + (diflon*diflat)
            ! average intensity
            clavint(INT(dat2d(x,y))) = clavint(INT(dat2d(x,y))) + pdat2d(x,y)
            ! peak intensity
            if(clpint(INT(dat2d(x,y)))<pdat2d(x,y))clpint(INT(dat2d(x,y))) = pdat2d(x,y)
          end do
        end do
        deallocate(dat2d,pdat2d)

      end do

      write(*,*)"======================================="
      write(*,*)"=== calculating (weighted) center of mass ..."
      write(*,*)"---------"

      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
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

        ! how many cells are in this time step?
        tp=0
        minclID=huge(minclID)
        maxclID=0
        maxarea=0
        do clID=1,globnIDs
          if(tsclID(clID)==tsID+1)then
            tp=tp+1
            if(clID<minclID)minclID=clID
            if(clID>maxclID)maxclID=clID
            if(clareagrd(clID)>maxarea)maxarea=clareagrd(clID)
          end if
        end do

        ! allocate array for coordinates for all gridpoints of cells in this tstep
        allocate(coords(minclID:maxclID,maxarea,2))
        allocate(coordsgrd(minclID:maxclID,maxarea,2))
        allocate(cellcnt(minclID:maxclID))
        allocate(ints(minclID:maxclID,maxarea))
        coords=-1
        coordsgrd=-1
        cellcnt=0
        ints=0

        ! now loop dat2d and calculate statistics
        do y=1,ny
          do x=1,nx
            if(dat2d(x,y)==outmissval)cycle
            ! assign coordinates for center of mass calculation later on
            ! save intensities to calculate weighted center of mass
            cellcnt(INT(dat2d(x,y)))=cellcnt(INT(dat2d(x,y)))+1
            coords(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),1)=xvals2d(x,y)
            coords(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),2)=yvals2d(x,y)
            coordsgrd(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),1)=x
            coordsgrd(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),2)=y
            ints(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))))=pdat2d(x,y)
          end do
        end do
        deallocate(dat2d,pdat2d)

        ! calculate center of masses
        if(periodic)then
          do clID=minclID,maxclID
            !first, find out if any grid point is at the boundary
            clbnd=.false.
            projx=.false.
            projy=.false.
            do i=1,cellcnt(clID)
              if(coords(clID,i,1)==minx .OR. coords(clID,i,1)==maxx .OR. &
                 & coords(clID,i,2)==miny .OR. coords(clID,i,2)==maxy)then
                clbnd=.true.
                exit
              end if
            end do
            do i=1,cellcnt(clID)
              if(coords(clID,i,1)==minx .OR. coords(clID,i,1)==maxx)projx=.true.
              if(coords(clID,i,2)==miny .OR. coords(clID,i,2)==maxy)projy=.true.
              if(projx .AND. projy)exit
            end do
            ! if that's the case count how many are closer to the upper boundary
            if(clbnd)then
              nrupbndx=0
              nrupbndy=0
              ! x direction
              if(projx)then
                do i=1,cellcnt(clID)
                  if(coords(clID,i,1)>(maxx/2))then
                    nrupbndx=nrupbndx+1
                  end if
                end do
                ! check at which boundary we want to project; x=1 or x=nx
                if((nrupbndx)>(clareagrd(clID)/2))then ! project at nx
                  if(verbose)write(*,*)"projecting cell",clID,"at nx with ",nrupbndx,"of",clareagrd(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,1)<(maxx/2))then
                      coords(clID,i,1)=coords(clID,i,1)+maxx
                    end if
                  end do                  
                else ! project at x=1
                  if(verbose)write(*,*)"projecting cell",clID,"at x=1 with ",clareagrd(clID)-nrupbndx,"of",clareagrd(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,1)>=(maxx/2))then
                      coords(clID,i,1)=minx-(maxx-coords(clID,i,1))
                    end if
                  end do                   
                end if
              end if
              ! y direction
              if(projy)then
                do i=1,cellcnt(clID)
                  if(coords(clID,i,2)>(maxy/2))then
                    nrupbndy=nrupbndy+1
                  end if
                end do
                ! check at which boundary we want to project; y=1 or y=ny
                if((nrupbndy)>(clareagrd(clID)/2))then ! project at ny
                  if(verbose)write(*,*)"projecting cell",clID,"at ny with ",nrupbndy,"of",clareagrd(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,2)<(maxy/2))then
                      coords(clID,i,2)=coords(clID,i,2)+maxy
                    end if
                  end do                  
                else ! project at y=1
                  if(verbose)write(*,*)"projecting cell",clID,"at y=1 with ",clareagrd(clID)-nrupbndy,"of",clareagrd(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,2)>=(maxy/2))then
                      coords(clID,i,2)=miny-(maxy-coords(clID,i,2))
                    end if
                  end do                   
                end if
              end if
              ! now calculate center of mass
              do i=1,cellcnt(clID)
                ! center of mass; on unit coordinate level
                clcmass(clID,1) = clcmass(clID,1) + coords(clID,i,1)
                clcmass(clID,2) = clcmass(clID,2) + coords(clID,i,2)
                wclcmass(clID,1) = wclcmass(clID,1) + coords(clID,i,1) * ints(clID,i)
                wclcmass(clID,2) = wclcmass(clID,2) + coords(clID,i,2) * ints(clID,i)
                ! center of mass; on grid point level
                clcmassgrd(clID,1) = clcmassgrd(clID,1) + coordsgrd(clID,i,1)
                clcmassgrd(clID,2) = clcmassgrd(clID,2) + coordsgrd(clID,i,2)
                wclcmassgrd(clID,1) = wclcmassgrd(clID,1) + coordsgrd(clID,i,1) * ints(clID,i)
                wclcmassgrd(clID,2) = wclcmassgrd(clID,2) + coordsgrd(clID,i,2) * ints(clID,i)
                wsum(clID) = wsum(clID) + ints(clID,i)
              end do              
            else ! if no gridpoint is at the boundaries, just do it the normal way
              do i=1,cellcnt(clID)
                ! center of mass; on unit coordinate level
                clcmass(clID,1) = clcmass(clID,1) + coords(clID,i,1)
                clcmass(clID,2) = clcmass(clID,2) + coords(clID,i,2)
                wclcmass(clID,1) = wclcmass(clID,1) + coords(clID,i,1) * ints(clID,i)
                wclcmass(clID,2) = wclcmass(clID,2) + coords(clID,i,2) * ints(clID,i)
                ! center of mass; on grid point level
                clcmassgrd(clID,1) = clcmassgrd(clID,1) + coordsgrd(clID,i,1)
                clcmassgrd(clID,2) = clcmassgrd(clID,2) + coordsgrd(clID,i,2)
                wclcmassgrd(clID,1) = wclcmassgrd(clID,1) + coordsgrd(clID,i,1) * ints(clID,i)
                wclcmassgrd(clID,2) = wclcmassgrd(clID,2) + coordsgrd(clID,i,2) * ints(clID,i)
                wsum(clID) = wsum(clID) + ints(clID,i)
              end do
            end if
          end do
        else ! no periodic boundaries
          do clID=minclID,maxclID
            do i=1,cellcnt(clID)
              ! center of mass; on unit coordinate level
              clcmass(clID,1) = clcmass(clID,1) + coords(clID,i,1)
              clcmass(clID,2) = clcmass(clID,2) + coords(clID,i,2)
              wclcmass(clID,1) = wclcmass(clID,1) + coords(clID,i,1) * ints(clID,i)
              wclcmass(clID,2) = wclcmass(clID,2) + coords(clID,i,2) * ints(clID,i)
              ! center of mass; on grid point level
              clcmassgrd(clID,1) = clcmassgrd(clID,1) + coordsgrd(clID,i,1)
              clcmassgrd(clID,2) = clcmassgrd(clID,2) + coordsgrd(clID,i,2)
              wclcmassgrd(clID,1) = wclcmassgrd(clID,1) + coordsgrd(clID,i,1) * ints(clID,i)
              wclcmassgrd(clID,2) = wclcmassgrd(clID,2) + coordsgrd(clID,i,2) * ints(clID,i)
              wsum(clID) = wsum(clID) + ints(clID,i)
            end do
          end do
        end if
        deallocate(coords,coordsgrd)
        deallocate(cellcnt)
        deallocate(ints)

      end do

      ! divide by area
      do i=1,globnIDs
        ! average intensity
        clavint(i)=clavint(i)/clareagrd(i)
        ! center of mass
        clcmass(i,1)=clcmass(i,1)/clareagrd(i)
        clcmass(i,2)=clcmass(i,2)/clareagrd(i)
        clcmassgrd(i,1)=clcmassgrd(i,1)/clareagrd(i)
        clcmassgrd(i,2)=clcmassgrd(i,2)/clareagrd(i)
        ! weighted center of mass
        wclcmass(i,1)=wclcmass(i,1)/wsum(i)
        wclcmass(i,2)=wclcmass(i,2)/wsum(i)
        wclcmassgrd(i,1)=wclcmassgrd(i,1)/wsum(i)
        wclcmassgrd(i,2)=wclcmassgrd(i,2)/wsum(i)
        ! bring projected center of mass back into the domain
        if(periodic)then
          ! in unit coordinates
          if(clcmass(i,1)>maxx)clcmass(i,1)=clcmass(i,1)-maxx
          if(clcmass(i,1)<minx)clcmass(i,1)=maxx-abs(clcmass(i,1))
          if(clcmass(i,2)>maxy)clcmass(i,2)=clcmass(i,2)-maxy
          if(clcmass(i,2)<miny)clcmass(i,2)=maxy-abs(clcmass(i,2))
          if(wclcmass(i,1)>maxx)wclcmass(i,1)=wclcmass(i,1)-maxx
          if(wclcmass(i,1)<minx)wclcmass(i,1)=maxx-abs(wclcmass(i,1))
          if(wclcmass(i,2)>maxy)wclcmass(i,2)=wclcmass(i,2)-maxy
          if(wclcmass(i,2)<miny)wclcmass(i,2)=maxy-abs(wclcmass(i,2))
        end if
      end do

      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,'(A,1i12)')  " AVerage cell area(gridpoints):",SUM(clareagrd)/globnIDs
      write(*,'(A,1f16.2)')" AVerage cell area(units)     :",SUM(clarea)/globnIDs
      write(*,'(A,1i12)')" Minimum cell area(gridpoints):",MINVAL(clareagrd)
      write(*,'(A,1f16.2)')" Minimum cell area(units)     :",MINVAL(clarea)
      write(*,'(A,1i12)')" Maximum cell area(gridpoints):",MAXVAL(clareagrd)
      write(*,'(A,1f16.2)')" Maximum cell area(units)     :",MAXVAL(clarea)
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
      write(1,*)"     clID    tsclID             clarea grd_clarea         clcmassX"//&
      &"         clcmassY        wclcmassX        wclcmassY     grd_clcmassX"//&
      &"     grd_clcmassY    grd_wclcmassX    grd_wclcmassY"//&
      &"            peakVal"//&
      &"              avVal  touchb      date(YYYYMMDD)        time(hhmmss)"
      do i=1,globnIDs
        write(1,'(2i10,1f19.2,1i11,8f17.3,2f19.10,1L8,2i20.6)')clIDs(i),tsclID(i),clarea(i),clareagrd(i),clcmass(i,:), &
          & wclcmass(i,:),clcmassgrd(i,:),wclcmassgrd(i,:),clpint(i),clavint(i), &
          & touchb(i),cldate(i),cltime(i)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"========= FINISHED STATISTICS ========="
      write(*,*)"======================================="

    end subroutine calccellstatistics

end module cellstatistics
