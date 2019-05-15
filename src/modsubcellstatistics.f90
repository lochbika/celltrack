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
module subcellstatistics

  use globvar

  implicit none

  contains
    subroutine calcsubcellstatistics

      use ncdfpars
      use globvar

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=stdfloattype), allocatable :: dat(:),pdat(:)          ! array for reading float from nc
      real(kind=stdfloattype), allocatable :: dat2d(:,:),pdat2d(:,:)  ! array for doing the clustering

      real(kind=stdfloattype), allocatable :: coords(:,:,:)     ! array for holding all cells coordinates
      real(kind=stdfloattype), allocatable :: ints(:,:)         ! array for holding all cells coordinates
      integer, allocatable :: cellcnt(:)             ! array for holding a counter for each cell

      logical :: clbnd ! true if a cell touches domain boundaries
      logical :: projx,projy ! shall we project points across a boundary
      integer :: nrupbndx,nrupbndy ! how many grid points are closer to the upper boundaries

      integer :: minclID,maxclID,maxarea

      real(kind=stdfloattype), allocatable :: wsum(:)

      write(*,*)"======================================="
      write(*,*)"========= SUBCELL STATISTICS =========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== collecting timesteps ..."
      write(*,*)"---------"

      !!!!!!!!!
      ! number of unique IDs is globnIDs
      ! loop time steps to find all IDs and their timesteps

      ! Ids range from 1 to globsubnIDs
      allocate(subclIDs(globsubnIDs))
      subclIDs=-1
      do i=1,globsubnIDs
        subclIDs(i)=i
      end do

      ! Open the cells file
      streamID2=streamOpenRead(suboutfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=getVarIDbyName(suboutfile,"subcellID")
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      allocate(subtsclID(globsubnIDs))
      allocate(subcldate(globsubnIDs))
      allocate(subcltime(globsubnIDs))
      subtsclID=-1
      tp=1
      ltp=1
      subcldate=-1
      subcltime=-1

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
          if(verbose)write(*,*)"NO subclusters in timestep:  ",tsID+1
          deallocate(dat)
          cycle
        end if

        ! assign time steps
        do i=1,nx*ny
          if(dat(i).ne.outmissval)then
            if(subtsclID(INT(dat(i)))==-1)then
              subtsclID(INT(dat(i)))=tsID+1
              subcldate(INT(dat(i)))=vdate(tsID+1)
              subcltime(INT(dat(i)))=vtime(tsID+1)
            end if
            if(verbose)write(*,*)"subcluster: ",subclIDs(INT(dat(i)))," is at timestep: ",subtsclID(INT(dat(i)))
            if(verbose)write(*,*)subcldate(INT(dat(i))),subcltime(INT(dat(i)))
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
      streamID2=streamOpenRead(suboutfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=getVarIDbyName(suboutfile,"subcellID")
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
      allocate(subclarea(globsubnIDs))
      allocate(subtouchb(globsubnIDs))
      allocate(subclcmass(globsubnIDs,2))
      allocate(subwclcmass(globsubnIDs,2))
      allocate(subclavint(globsubnIDs))
      allocate(subclpint(globsubnIDs))
      allocate(wsum(globsubnIDs))
      subclarea=0
      subtouchb=.false.
      subclavint=0
      subclpint=0
      wsum=0
      subclcmass=0
      subwclcmass=0

      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if

        ! if this time step has no values: set touchb to true for all cells in previous ts
        ! do this also for cells in next time step
        ! because it could be that this interrupts tracks in the middle of their life time
        if(tsALLna(tsID+1))then
          do i=1,globsubnIDs
            if(subtsclID(i)==tsID .OR. subtsclID(i)==tsID+2)subtouchb(i)=.true.
          end do
        end if

        ! cycle if there are no cells in this time step
        if(.NOT.ANY(subtsclID==tsID+1))cycle

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
              if(pdat2d(x-1,y)==inmissval)subtouchb(INT(dat2d(x,y))) = .true.
            end if
            if(y.ne.1)then
              if(pdat2d(x,y-1)==inmissval)subtouchb(INT(dat2d(x,y))) = .true.
            end if
            if(x.ne.nx)then
              if(pdat2d(x+1,y)==inmissval)subtouchb(INT(dat2d(x,y))) = .true.
            end if
            if(y.ne.ny)then
              if(pdat2d(x,y+1)==inmissval)subtouchb(INT(dat2d(x,y))) = .true.
            end if
            ! this cell touches time or space boundaries?
            if(periodic)then
              if(tsID==0 .OR. tsID==ntp-1)subtouchb(INT(dat2d(x,y))) = .true.
            else
              if(y==1 .OR. x==1 .OR. x==nx .OR. y==ny .OR. tsID==0 .OR. tsID==ntp-1)subtouchb(INT(dat2d(x,y))) = .true.
            end if
            ! cell area
            subclarea(INT(dat2d(x,y))) = subclarea(INT(dat2d(x,y))) + 1
            ! average intensity
            subclavint(INT(dat2d(x,y))) = subclavint(INT(dat2d(x,y))) + pdat2d(x,y)
            ! peak intensity
            if(subclpint(INT(dat2d(x,y)))<pdat2d(x,y))subclpint(INT(dat2d(x,y))) = pdat2d(x,y)
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
        if(.NOT.ANY(subtsclID==tsID+1))cycle

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
        do clID=1,globsubnIDs
          if(subtsclID(clID)==tsID+1)then
            tp=tp+1
            if(clID<minclID)minclID=clID
            if(clID>maxclID)maxclID=clID
            if(subclarea(clID)>maxarea)maxarea=subclarea(clID)
          end if
        end do

        ! allocate array for coordinates for all gridpoints of cells in this tstep
        allocate(coords(minclID:maxclID,maxarea,2))
        allocate(cellcnt(minclID:maxclID))
        allocate(ints(minclID:maxclID,maxarea))
        coords=-1
        cellcnt=0
        ints=0

        ! now loop dat2d and calculate statistics
        do y=1,ny
          do x=1,nx
            if(dat2d(x,y)==outmissval)cycle
            ! assign coordinates for center of mass calculation later on
            ! save intensities to calculate weighted center of mass
            cellcnt(INT(dat2d(x,y)))=cellcnt(INT(dat2d(x,y)))+1
            coords(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),1)=x
            coords(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),2)=y
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
              if(coords(clID,i,1)==1 .OR. coords(clID,i,1)==nx .OR. &
                 & coords(clID,i,2)==1 .OR. coords(clID,i,2)==ny)then
                clbnd=.true.
                exit
              end if
            end do
            do i=1,cellcnt(clID)
              if(coords(clID,i,1)==1 .OR. coords(clID,i,1)==nx)projx=.true.
              if(coords(clID,i,2)==1 .OR. coords(clID,i,2)==ny)projy=.true.
              if(projx .AND. projy)exit
            end do
            ! if that's the case count how many are closer to the upper boundary
            if(clbnd)then
              nrupbndx=0
              nrupbndy=0
              ! x direction
              if(projx)then
                do i=1,cellcnt(clID)
                  if(coords(clID,i,1)>(nx/2))then
                    nrupbndx=nrupbndx+1
                  end if
                end do
                ! check at which boundary we want to project; x=1 or x=nx
                if((nrupbndx)>(subclarea(clID)/2))then ! project at nx
                  if(verbose)write(*,*)"projecting cell",clID,"at nx with ",nrupbndx,"of",subclarea(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,1)<(nx/2))then
                      coords(clID,i,1)=coords(clID,i,1)+nx
                    end if
                  end do                  
                else ! project at x=1
                  if(verbose)write(*,*)"projecting cell",clID,"at x=1 with ",subclarea(clID)-nrupbndx,"of",subclarea(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,1)>=(nx/2))then
                      coords(clID,i,1)=1-(nx-coords(clID,i,1))
                    end if
                  end do                   
                end if
              end if
              ! y direction
              if(projy)then
                do i=1,cellcnt(clID)
                  if(coords(clID,i,2)>(ny/2))then
                    nrupbndy=nrupbndy+1
                  end if
                end do
                ! check at which boundary we want to project; y=1 or y=ny
                if((nrupbndy)>(subclarea(clID)/2))then ! project at ny
                  if(verbose)write(*,*)"projecting cell",clID,"at ny with ",nrupbndy,"of",subclarea(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,2)<(ny/2))then
                      coords(clID,i,2)=coords(clID,i,2)+ny
                    end if
                  end do                  
                else ! project at y=1
                  if(verbose)write(*,*)"projecting cell",clID,"at y=1 with ",subclarea(clID)-nrupbndy,"of",subclarea(clID),"points"
                  do i=1,cellcnt(clID)
                    if(coords(clID,i,2)>=(ny/2))then
                      coords(clID,i,2)=1-(ny-coords(clID,i,2))
                    end if
                  end do                   
                end if
              end if
              ! now calculate center of mass
              do i=1,cellcnt(clID)
                ! center of masses
                subclcmass(clID,1) = subclcmass(clID,1) + coords(clID,i,1)
                subclcmass(clID,2) = subclcmass(clID,2) + coords(clID,i,2)
                subwclcmass(clID,1) = subwclcmass(clID,1) + coords(clID,i,1) * ints(clID,i)
                subwclcmass(clID,2) = subwclcmass(clID,2) + coords(clID,i,2) * ints(clID,i)
                wsum(clID) = wsum(clID) + ints(clID,i)
              end do              
            else ! if no gridpoint is at the boundaries, just do it the normal way
              do i=1,cellcnt(clID)
                ! normal center of masses
                subclcmass(clID,1) = subclcmass(clID,1) + coords(clID,i,1)
                subclcmass(clID,2) = subclcmass(clID,2) + coords(clID,i,2)
                subwclcmass(clID,1) = subwclcmass(clID,1) + coords(clID,i,1) * ints(clID,i)
                subwclcmass(clID,2) = subwclcmass(clID,2) + coords(clID,i,2) * ints(clID,i)
                wsum(clID) = wsum(clID) + ints(clID,i)
              end do
            end if
          end do
        else ! no periodic boundaries
          do clID=minclID,maxclID
            do i=1,cellcnt(clID)
              ! normal center of masses
              subclcmass(clID,1) = subclcmass(clID,1) + coords(clID,i,1)
              subclcmass(clID,2) = subclcmass(clID,2) + coords(clID,i,2)
              subwclcmass(clID,1) = subwclcmass(clID,1) + coords(clID,i,1) * ints(clID,i)
              subwclcmass(clID,2) = subwclcmass(clID,2) + coords(clID,i,2) * ints(clID,i)
              wsum(clID) = wsum(clID) + ints(clID,i)
            end do
          end do
        end if
        deallocate(coords)
        deallocate(cellcnt)
        deallocate(ints)

      end do

      ! divide by area
      do i=1,globsubnIDs
        ! average intensity
        subclavint(i)=subclavint(i)/subclarea(i)
        ! center of mass
        subclcmass(i,1)=subclcmass(i,1)/subclarea(i)
        subclcmass(i,2)=subclcmass(i,2)/subclarea(i)
        ! weighted center of mass
        subwclcmass(i,1)=subwclcmass(i,1)/wsum(i)
        subwclcmass(i,2)=subwclcmass(i,2)/wsum(i)
        ! bring projected center of mass back into the domain
        if(subclcmass(i,1)>nx)subclcmass(i,1)=subclcmass(i,1)-nx
        if(subclcmass(i,1)<1)subclcmass(i,1)=nx-abs(subclcmass(i,1))
        if(subclcmass(i,2)>ny)subclcmass(i,2)=subclcmass(i,2)-ny
        if(subclcmass(i,2)<1)subclcmass(i,2)=ny-abs(subclcmass(i,2))
        if(subwclcmass(i,1)>nx)subwclcmass(i,1)=subwclcmass(i,1)-nx
        if(subwclcmass(i,1)<1)subwclcmass(i,1)=nx-abs(subwclcmass(i,1))
        if(subwclcmass(i,2)>ny)subwclcmass(i,2)=subwclcmass(i,2)-ny
        if(subwclcmass(i,2)<1)subwclcmass(i,2)=ny-abs(subwclcmass(i,2))
      end do

      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,'(A,1i12)')" AVerage subcell area(gridpoints):",SUM(subclarea)/globsubnIDs
      write(*,'(A,1i12)')" Minimum subcell area(gridpoints):",MINVAL(subclarea)
      write(*,'(A,1i12)')" Maximum subcell area(gridpoints):",MAXVAL(subclarea)
      write(*,'(A,1f12.6)')" TOtal average value           :",SUM(subclavint)/globsubnIDs
      write(*,'(A,1f12.6)')" Mininmum avergae value        :",MINVAL(subclavint)
      write(*,'(A,1f12.6)')" Maxinmum avergae value        :",MAXVAL(subclavint)
      write(*,'(A,1f12.6)')" AVerage peak value            :",SUM(subclpint)/globsubnIDs
      write(*,'(A,1f12.6)')" Mininmum peak value           :",MINVAL(subclpint)
      write(*,'(A,1f12.6)')" Maxinmum peak value           :",MAXVAL(subclpint)
      write(*,*)"---------"

      write(*,*)"======================================="
      write(*,*)"=== writing stats to file subcell_stats.txt ..."
      write(*,*)"---------"

      ! write the stats to file
      open(unit=1,file="subcell_stats.txt",action="write",status="replace")
      write(1,*)"   subclID  subtsclID  subclarea     subclcmassX"//&
      &"     subclcmassY    subwclcmassX    subwclcmassY          subpeakVal"//&
      &"            subavVal subtouchb       date(YYYYMMDD)         time(hhmmss)"
      do i=1,globsubnIDs
        write(1,'(3i11,4f16.6,2f20.12,1L10,2i21.6)')subclIDs(i),subtsclID(i),subclarea(i),subclcmass(i,:), &
          & subwclcmass(i,:),subclpint(i),subclavint(i),subtouchb(i),subcldate(i),subcltime(i)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"========= FINISHED STATISTICS ========="
      write(*,*)"======================================="

    end subroutine calcsubcellstatistics

end module subcellstatistics
