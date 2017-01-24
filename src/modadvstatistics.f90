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
      real(kind=8) :: probs(22)

      write(*,*)"======================================="
      write(*,*)"===== ADDITIONAL CELL STATISTICS ======"
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== calculating percentiles for all cells ..."
      write(*,*)"---------"

      ! prepare probabilities for quantiles
      probs=(/0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.975,0.99,0.999/)

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
        do k=1,nx*ny
          if(dat(k)==outmissval)cycle
          cellvalues( INT(dat(k)) , cellcounter(INT(dat(k))) )=pdat(k)
          cellcounter(INT(dat(k)))=cellcounter(INT(dat(k)))+1
        end do
        deallocate(dat,pdat)
      end do
      deallocate(cellcounter)

      CALL streamClose(streamID1)
      CALL streamClose(streamID2)
      
      ! now sort each cells values in ascending order
      do i=1,globnIDs
        CALL QsortC(cellvalues(i,1:clarea(i)))
      end do
      
      ! now calculate cell value percentiles
      allocate(cellperc(globnIDs,24))
      cellperc=0.D0
      do i=1,globnIDs
        cellperc(i,1)=MINVAL(cellvalues(i,1:clarea(i)))
        cellperc(i,24)=MAXVAL(cellvalues(i,1:clarea(i)))
        do p=1,22
          CALL quantile(cellvalues(i,1:clarea(i)),probs(p),cellperc(i,p+1))
        end do
      end do      

      write(*,*)"======================================="
      write(*,*)"=== writing stats to file cell_advstats_percentiles.txt ..."
      write(*,*)"---------"

      ! write the stats to file
      open(unit=1,file="cell_advstats_percentiles.txt",action="write",status="replace")
      write(1,*)"     clID            MIN           0.05           0.10           0.15"//&
      &"           0.20           0.25           0.30"//&
      &"           0.35           0.40           0.45"//&
      &"           0.50           0.55           0.60"//&
      &"           0.65           0.70           0.75"//&
      &"           0.80           0.85           0.90"//&
      &"           0.95          0.975           0.99"//&
      &"          0.999            MAX"
      do i=1,globnIDs
        write(1,'(1i10,24f15.6)')clIDs(i),cellperc(i,:)
      end do
      close(unit=1)

      deallocate(cellperc,cellvalues)

      write(*,*)"======================================="
      write(*,*)"=== FINISHED ADDITIONAL STATISTICS ===="
      write(*,*)"======================================="

    end subroutine calccellpercentiles
    
    subroutine quantile(x,prob,qu)
      ! taken from: http://fortranwiki.org/fortran/show/Quartiles
      implicit none
      real(kind=8), intent(in)  :: x(:)
      real(kind=8), intent(in)  :: prob
      real(kind=8), intent(out) :: qu
      
      integer :: n,ib
      real(kind=8) :: tol,a,b,c,diff
      
      n=SIZE(x)
      tol=1.e-8

      a=n*prob
      CALL getgp(a,b,c)
      
      ib=int(c)
     
      diff=b-0.0D0
      if(diff<=tol) then
        qu=(x(ib+1)+x(ib))/2.0d0
      else
        qu=x(ib+1)
      end if
    end subroutine quantile
    
    subroutine getgp(a,b,c)
      ! taken from: http://fortranwiki.org/fortran/show/Quartiles
      implicit none
      ! Subroutine to that returns the Right hand and Left hand side digits of a decimal number
      real(kind=8) :: a,b,c
      
      b=mod(a,1.0d0)
      c=a-b
    end subroutine getgp

end module advstats
