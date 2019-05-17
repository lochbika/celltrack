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
module trackstatistics

  use globvar

  implicit none

  contains
    subroutine calctrackstatistics

      use globvar

      implicit none

      character(len=stdclen) :: ttrack
      real(kind=stdfloattype) :: pint,avint
      integer :: pinttime

      write(*,*)"======================================="
      write(*,*)"========== TRACK STATISTICS ==========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Reading tracks_all.txt ..."
      write(*,*)"---------"

      ! scan tracks_all_stats.txt
      CALL scan_trackfile("tracks_all.txt",ntracks,maxtrlen)

      ! open the tracks file
      open(unit=1,file="tracks_all.txt",action="read",status="old")

      !now loop all lines in tracks_all.txt
      ctrack=0
      ! save all tracks in this buffer
      allocate(alltracks(ntracks,maxtrlen),trtypes(ntracks))
      allocate(nobounds(ntracks),trpint(ntracks),trpinttime(ntracks),travint(ntracks),trdur(ntracks))
      alltracks=-1

      do
        ! read the line ;)
        read(1,'(A)',IOSTAT=riostat)ttrack
        if(riostat.ne.0)exit
        !ttrack=trim(ttrack)
        !check if a new track begins
        if(ttrack(1:3)=="###")then
          read(1,*)
          ! buffer counter
          tp=1
          read(ttrack(5:16),'(1i12)')ctrack
          read(ttrack(20:20),'(1L1)')nobounds(ctrack)
          read(ttrack(21:24),'(1i4)')trtypes(ctrack)
          if(verbose)then
            if(MOD(ctrack,outstep)==0 .OR. ctrack==1 .OR. ctrack==ntracks)then
              write(*,*)"Reading track with number/ID ",ctrack," of ",ntracks
            end if
          end if
        else
          !now read cell ID, time step and area
          read(ttrack,'(1i12)')alltracks(ctrack,tp)
          tp=tp+1
        end if
      end do

      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"=== Calculating track statistics ..."
      write(*,*)"---------"

      ! open the track stats and summary files
      open(unit=3,file="tracks_all_summary.txt",action="write",status="replace")
      open(unit=4,file="tracks_all_stats.txt",action="write",status="replace")
      open(unit=2,file="tracks_clean_summary.txt",action="write",status="replace")
      open(unit=5,file="tracks_clean_stats.txt",action="write",status="replace")
      ! write header
      write(2,*)"    trackID      trType            peakVal    pValtime              avVal start   dur"
      write(3,*)"    trackID      trType            peakVal    pValtime              avVal start   dur"

      ! loop all tracks
      do i=1,ntracks
        ! status
        if(MOD(i,outstep)==0 .OR. i==1 .OR. i==ntracks)then
          write(*,*)"Processing track with number/ID ",i," of ",ntracks
        end if

        ! write stats
        write(4,'(1a4,1i12,1L4,1i4)')"### ",i,nobounds(i),trtypes(i)
        write(4,*)"     clID    tsclID             clarea grd_clarea         clcmassX"//&
        &"         clcmassY        wclcmassX        wclcmassY     grd_clcmassX"//&
        &"     grd_clcmassY     grd_clcmassX     grd_clcmassY"//&
        &"            peakVal"//&
        &"              avVal  touchb      date(YYYYMMDD)        time(hhmmss)"
        do k=1,maxtrlen
          if(alltracks(i,k)==-1)exit
          write(4,'(2i10,1f19.2,1i11,8f17.3,2f19.10,1L8,2i20.6)')clIDs(alltracks(i,k)),tsclID(alltracks(i,k)), &
            & clarea(alltracks(i,k)),clareagrd(alltracks(i,k)),clcmass(alltracks(i,k),:), &
            & wclcmass(alltracks(i,k),:),clcmassgrd(alltracks(i,k),:),wclcmassgrd(alltracks(i,k),:),&
            & clpint(alltracks(i,k)),clavint(alltracks(i,k)),touchb(alltracks(i,k))
        end do
        ! write stats for clean tracks
        if(nobounds(i) .AND. trtypes(i)==9)then
          write(5,'(1a4,1i12,1L4,1i4)')"### ",i,nobounds(i),trtypes(i)
          write(5,*)"     clID    tsclID             clarea grd_clarea         clcmassX"//&
        &"         clcmassY        wclcmassX        wclcmassY     grd_clcmassX"//&
        &"     grd_clcmassY     grd_clcmassX     grd_clcmassY"//&
        &"            peakVal"//&
        &"              avVal  touchb      date(YYYYMMDD)        time(hhmmss)"
          do k=1,maxtrlen
            if(alltracks(i,k)==-1)exit
            write(5,'(2i10,1f19.2,1i11,8f17.3,2f19.10,1L8,2i20.6)')clIDs(alltracks(i,k)),tsclID(alltracks(i,k)), &
            & clarea(alltracks(i,k)),clareagrd(alltracks(i,k)),clcmass(alltracks(i,k),:), &
            & wclcmass(alltracks(i,k),:),clcmassgrd(alltracks(i,k),:),wclcmassgrd(alltracks(i,k),:),&
            & clpint(alltracks(i,k)),clavint(alltracks(i,k)),touchb(alltracks(i,k))
          end do
        end if

        ! calc summary values
        avint=0
        pint=0
        tp=0
        do k=1,maxtrlen
          if(alltracks(i,k)==-1)exit
          if(clpint(alltracks(i,k))>pint)then
            pint=clpint(alltracks(i,k))
            pinttime=tsclID(alltracks(i,k))
          end if
          avint=avint+clavint(alltracks(i,k))
          tp=tp+1
        end do
        travint(i)=avint/tp
        trpint(i)=pint
        trpinttime(i)=pinttime
        trdur(i)=tsclID(alltracks(i,tp))-tsclID(alltracks(i,1))+1

        ! write to file
        write(3,'(2i12,1f19.12,1i12,1f19.12,2i6)')i,trtypes(i),trpint(i),trpinttime(i), &
          & travint(i),tsclID(alltracks(i,1)),trdur(i)

        ! write to file for clean tracks
        if(nobounds(i) .AND. trtypes(i)==9)then
          write(2,'(2i12,1f19.12,1i12,1f19.12,2i6)')i,trtypes(i),trpint(i),trpinttime(i), &
            & travint(i),tsclID(alltracks(i,1)),trdur(i)
        end if

      end do

      close(unit=1)
      close(unit=3)
      close(unit=4)
      close(unit=2)
      close(unit=5)
      
      write(*,*)"======================================="
      write(*,*)"====== FINISHED TRACK STATISTICS ======"
      write(*,*)"======================================="

    end subroutine calctrackstatistics

    subroutine writetracks

      use ncdfpars
      use globvar

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=stdfloattype), allocatable :: dat(:)      ! array for reading float from nc

      integer, allocatable :: tnums(:)
      
      write(*,*)"======================================="
      write(*,*)"========= WRITE TRACKS TO FILE ========"
      write(*,*)"======================================="
      write(*,*)"=== ... to tracks.nc ..."
      write(*,*)"---------"

      ! Open the cells file
      streamID2=streamOpenRead(outfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if

      ! Set the IDs
      varID2=0
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)

      !! open new nc file for results
      ! define grid
      gridID1=gridDuplicate(gridID2)
      zaxisID1=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID1, level)
      ! define variables
      vlistID1=vlistCreate()
      varID1=vlistDefVar(vlistID1,gridID1,zaxisID1,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID1,varID1,"trackID")
      CALL vlistDefVarLongname(vlistID1,varID1,"unique ID of each track")
      CALL vlistDefVarUnits(vlistID1,varID1,"-")
      CALL vlistDefVarMissval(vlistID1,varID1,outmissval)
      CALL vlistDefVarDatatype(vlistID1,varID1,CDI_DATATYPE_INT32)
      ! copy time axis from input
      taxisID1=vlistInqTaxis(vlistID2)
      call vlistDefTaxis(vlistID1,taxisID1)
      ! Open the dataset for writing
      streamID1=streamOpenWrite("tracks.nc",CDI_FILETYPE_NC4)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if
      ! set netCDF4 compression
      CALL streamDefCompType(streamID1,CDI_COMPRESS_ZIP)
      CALL streamDefCompLevel(streamID1, 6)
      ! Assign variables to dataset
      call streamDefVList(streamID1,vlistID1)

      ! allocate data arrays for storage
      allocate(dat(nx*ny))

      do tsID=0,ntp-1
        ! status to stdout
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1
        end if
      
        ! Set time step for input and output
        status=streamInqTimestep(streamID2,tsID)
        status=streamDefTimestep(streamID1,tsID)

        ! Read time step from input
        call streamReadVar(streamID2,varID2,dat,nmiss2)

        ! check which tracks are in this time step
        allocate(tnums(5000))
        tnums=-1
        tp=1
        do i=1,ntracks
          do j=1,maxtrlen
            if(alltracks(i,j)==-1)exit
            if(tsclID(alltracks(i,j))==tsID+1)then
              tnums(tp)=i
              tp=tp+1
              exit
            end if
          end do
        end do

        ! now loop dat and replace clIDs with trIDs
        do i=1,nx*ny
          if(dat(i)==outmissval)cycle
          do j=1,5000
            if(tnums(j)==-1)exit
            do k=1,maxtrlen
              if(alltracks(tnums(j),k)==-1)exit
              if(INT(dat(i))==alltracks(tnums(j),k))dat(i)=tnums(j)
            end do
          end do
        end do
        deallocate(tnums)

        ! write this time step to tracks.nc
        call streamWriteVar(streamID1,varID1,dat,nmiss1)

      end do

      deallocate(dat)
      CALL gridDestroy(gridID1)
      CALL vlistDestroy(vlistID1)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"======= FINISHED WRITING TRACKS ======="
      write(*,*)"======================================="

    end subroutine writetracks

end module trackstatistics
