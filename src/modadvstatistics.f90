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
module advstats

  use globvar
  use advstatsdata
  use quicksort

  implicit none

  contains
    subroutine calccellpercentiles

      use ncdfpars
      use globvar
      use advstatsdata

      implicit none

      include 'cdi.inc'

      ! local data arrays
      real(kind=8), allocatable :: dat(:),pdat(:)     ! array for reading float from nc
      real(kind=8), allocatable :: cellvalues(:,:)    ! holds all values of each cell
      integer, allocatable :: cellcounter(:)          ! current position for each cells grid points

      write(*,*)"======================================="
      write(*,*)"===== ADDITIONAL CELL STATISTICS ======"
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== calculating percentiles for all cells ..."
      write(*,*)"---------"

      !!!!!!!!!
      ! gather all cell values and calculate percentiles
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
      allocate(cellvalues(globnIDs,MAXVAL(clarea)))
      allocate(cellcounter(globnIDs))
      cellvalues=-1
      cellcounter=1

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

        ! now loop dat and gather each cells values
        do y=1,ny
          do x=1,nx
            if(dat(x,y)==outmissval)cycle
            cellvalues(INT(dat(x,y)),cellcounter(INT(dat(x,y))))=pdat(x,y)
          end do
        end do
        deallocate(dat,pdat)
      end do

      CALL streamClose(streamID1)
      CALL streamClose(streamID2)
      
      ! now sort each cells values in ascending order
      do i=1,globnIDs
        CALL QsortC(cellvalues(i,:))
      end do

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
      write(*,*)"=== FINISHED ADDITIONAL STATISTICS ===="
      write(*,*)"======================================="

    end subroutine calccellpercentiles

end module advstats
