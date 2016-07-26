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

      real(kind=8) :: wsum
      real(kind=8), allocatable :: clgrdc(:,:),wclgrdc(:,:)

      logical, allocatable :: mask1d(:)
      
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
      streamID1=streamOpenRead(outfile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
      varID1=0
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      allocate(tsclID(globnIDs))
      tsclID=-1
      tp=1
      ltp=1

      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if

        ! Allocate arrays for data storage
        allocate(dat(nx*ny))

        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)

        ! Read time step from input
        call streamReadVar(streamID1,varID1,dat,nmiss2)

        !cycle if all values are -999
        if(ALL(dat==-999.D0))then
          if(verbose)write(*,*)"NO clusters in timestep:  ",tsID+1
          deallocate(dat)
          cycle
        end if

        ! assign time steps
        do i=1,nx*ny
          if(dat(i).ne.-999.D0)then
            tsclID(INT(dat(i)))=tsID+1
            if(verbose)write(*,*)"cluster: ",clIDs(INT(dat(i)))," is at timestep: ",tsclID(INT(dat(i)))
          end if
        end do
        deallocate(dat)
      end do

      CALL streamClose(streamID1)

      write(*,*)"======================================="
      write(*,*)"=== calculating area, (weighted) center of mass, peak and average values ..."
      write(*,*)"---------"

      !!!!!!!!!
      ! find center of mass and area of clusters
      CALL datainfo(outfile)

      ! Open the cells file
      streamID1=streamOpenRead(outfile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
      varID1=0
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

      ! Open the original data file
      streamID2=streamOpenRead(ifile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      varID2=ivar
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      allocate(clarea(globnIDs))
      allocate(wclarea(globnIDs))
      allocate(touchb(globnIDs))
      touchb=.false.
      allocate(clcmass(globnIDs,2))
      allocate(wclcmass(globnIDs,2))
      allocate(clavint(globnIDs))
      allocate(clpint(globnIDs))
      clavint=0
      clpint=0


      do i=1,globnIDs
        if(MOD(i,outstep*5)==0 .OR. i==1 .OR. i==globnIDs)then
          write(*,*)"Processing cell: ",i,"/",globnIDs,"..."
        end if

        ! Allocate arrays for data storage
        allocate(dat(nx*ny))
        allocate(pdat(nx*ny))
        ! Allocate array for masking
        allocate(mask1d(nx*ny))
        mask1d=.false.

        ! set timestep according to cluster
        tsID=tsclID(i)-1

        ! get cluster ID
        clID=clIDs(i)

        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)
        status=streamInqTimestep(streamID2,tsID)

        ! Read time step from input
        call streamReadVar(streamID1,varID1,dat,nmiss2)
        call streamReadVarSlice(streamID2,varID2,levelID,pdat,nmiss2)

        ! mask gridpoints with this cluster
        do k=1,nx*ny
          if(dat(k)==clID)mask1d(k)=.true.
        end do

        ! now calc area of this cluster
        clarea(i)=SUM(dat,MASK=mask1d)/clID
        if(verbose)write(*,*)"cluster: ",clID," has area: ",clarea(i)

        !reashape to 2d for center of mass calculation
        allocate(dat2d(nx,ny),pdat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        CALL reshapeT2d(pdat,nx,ny,pdat2d)
        deallocate(dat,pdat)

        ! allocate coordinates array
        allocate(clgrdc(clarea(i),2))
        allocate(wclgrdc(clarea(i),2))

        ! now loop dat2d and check if gridpoints belong to cluster
        tp=1
        k=-1
        wsum=0
        do y=1,ny
          do x=1,nx
            if(dat2d(x,y)==clID)then
              if(y==1 .OR. x==1 .OR. x==nx .OR. y==ny)touchb(i)=.true.
              clgrdc(tp,1)=x
              clgrdc(tp,2)=y
              wclgrdc(tp,1)=x*pdat2d(x,y)
              wclgrdc(tp,2)=y*pdat2d(x,y)
              wsum=wsum+pdat2d(x,y)
              if(clpint(i)<pdat2d(x,y))clpint(i)=pdat2d(x,y)
              clavint(i)=clavint(i)+pdat2d(x,y)
              tp=tp+1
              k=0
            else
              if(k>-1)k=k+1
            end if
          end do
          if(k>2*nx)exit
        end do
        clavint(i)=clavint(i)/clarea(i)
        clcmass(i,1)=SUM(clgrdc(:,1))/clarea(i)
        clcmass(i,2)=SUM(clgrdc(:,2))/clarea(i)
        wclcmass(i,1)=SUM(wclgrdc(:,1))/wsum
        wclcmass(i,2)=SUM(wclgrdc(:,2))/wsum
        if(verbose)write(*,*)"cluster: ",clID," has center of mass at: ",clcmass(i,:)

        deallocate(mask1d,clgrdc,wclgrdc,dat2d,pdat2d)

      end do

      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,'(A,1i12)')"AVerage cell area(gridpoints):",SUM(clarea)/globnIDs
      write(*,'(A,1i12)')"Minimum cell area(gridpoints):",MINVAL(clarea)
      write(*,'(A,1i12)')"Maximum cell area(gridpoints):",MAXVAL(clarea)
      write(*,'(A,1f12.6)')"TOtal average value          :",SUM(clavint)/globnIDs
      write(*,'(A,1f12.6)')"Mininmum avergae value       :",MINVAL(clavint)
      write(*,'(A,1f12.6)')"Maxinmum avergae value       :",MAXVAL(clavint)
      write(*,'(A,1f12.6)')"AVerage peak value           :",SUM(clpint)/globnIDs
      write(*,'(A,1f12.6)')"Mininmum peak value          :",MINVAL(clpint)
      write(*,'(A,1f12.6)')"Maxinmum peak value          :",MAXVAL(clpint)
      write(*,'(A,1f12.6)')"---------"

      write(*,*)"---------"
      write(*,*)"======================================="
      write(*,*)"=== writing stats to file cell_stats.txt ..."
      write(*,*)"---------"

      ! write the stats to file
      open(unit=1,file="cell_stats.txt",action="write",status="replace")
      write(1,*)"     clID    tsclID    clarea       clcmassX"//&
      &"       clcmassY      wclcmassX      wclcmassY            peakVal              avVal  touchb"
      do i=1,globnIDs
        write(1,'(3i10,4f15.6,2f19.12,1L8)')clIDs(i),tsclID(i),clarea(i),clcmass(i,:), &
          & wclcmass(i,:),clpint(i),clavint(i),touchb(i)
      end do
      close(unit=1)
      
      write(*,*)"======================================="
      write(*,*)"========= FINISHED STATISTICS ========="
      write(*,*)"======================================="

    end subroutine calccellstatistics

end module cellstatistics
