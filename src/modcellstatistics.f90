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
      allocate(clIDs(globnIDs))
      clIDs=-1
      do i=1,globnIDs
        clIDs(i)=i
      end do

      CALL datainfo(outfile)

      ! Open the cells file
      streamID2=streamOpenRead(outfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=0
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
      write(*,*)"=== calculating area, (weighted) center of mass, peak and average values ..."
      write(*,*)"---------"

      !!!!!!!!!
      ! find center of mass and area of clusters
      CALL datainfo(outfile)

      ! Open the cells file
      streamID2=streamOpenRead(outfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=0
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
      varID1=ivar
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      ! allocate and initialize arrays for cell statistics
      allocate(clarea(globnIDs))
      allocate(touchb(globnIDs))
      allocate(clcmass(globnIDs,2))
      allocate(wclcmass(globnIDs,2))
      allocate(clavint(globnIDs))
      allocate(clpint(globnIDs))
      allocate(wsum(globnIDs))
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

    end subroutine calccellstatistics

end module cellstatistics
