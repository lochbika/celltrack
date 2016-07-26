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
module metatrackstatistics

  use globvar

  implicit none

  contains
    subroutine calcmetatrackstatistics

      use globvar

      implicit none

      character(len=800) :: ttrack
      integer :: mints,maxts

      write(*,*)"======================================="
      write(*,*)"======== META TRACK STATISTICS ========"
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Reading meta_all.txt ..."
      write(*,*)"---------"

      ! scan meta_all.txt
      CALL scan_trackfile("meta_all.txt",nmeta,maxmetalen)

      ! open the tracks file
      open(unit=1,file="meta_all.txt",action="read",status="old")

      !now loop all lines in tracks_all.txt
      cmeta=0
      ! save all tracks in this buffer
      allocate(allmeta(nmeta,maxmetalen),mnobounds(nmeta))
      allmeta=-1

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
          read(ttrack(5:16),'(1i12)')cmeta
          read(ttrack(20:20),'(1L1)')mnobounds(cmeta)
          if(verbose)then
            if(MOD(cmeta,outstep)==0 .OR. cmeta==1 .OR. cmeta==nmeta)then
              write(*,*)"Reading meta track with number/ID ",cmeta," of ",nmeta
            end if
          end if
        else
          !now read cell ID, time step and area
          read(ttrack,'(1i12)')allmeta(cmeta,tp)
          tp=tp+1
        end if
      end do

      close(unit=1)

      ! save to which meta track a cell belongs
      allocate(clmeta(globnIDs))
      clmeta=-1
      do i=1,nmeta
        do k=1,maxmetalen
          if(allmeta(i,k)==-1)exit
          do j=1,maxtrlen
            if(alltracks(allmeta(i,k),j)==-1)exit
            clmeta(alltracks(allmeta(i,k),j))=i
          end do
        end do
      end do

      write(*,*)"======================================="
      write(*,*)"=== meta tracks stats meta_stats.txt ..."
      write(*,*)"---------"

      open(unit=1,file="meta_stats.txt",action="write",status="replace")
      do i=1,nmeta
        ! write header meta_stats.txt
        write(1,'(1a4,1i12,1L4,1i4)')"### ",i,mnobounds(i)
        write(1,*)"    trackID      trType            peakVal    pValtime              avVal start   dur"
        do k=1,maxmetalen
          if(allmeta(i,k)==-1)exit
          write(1,'(2i12,1f19.12,1i12,1f19.12,2i6)')allmeta(i,k),trtypes(allmeta(i,k)), &
          & trpint(allmeta(i,k)),trpinttime(allmeta(i,k)), &
          & travint(allmeta(i,k)),tsclID(alltracks(i,1)),trdur(allmeta(i,k))
        end do
      end do

      write(*,*)"======================================="
      write(*,*)"=== meta tracks summary meta_summary.txt ..."
      write(*,*)"---------"

      ! loop meta tracks
      ! 10,12,17,18,20,33,34,36
      allocate(nkinds(nmeta,8))
      nkinds=0
      do i=1,nmeta
        ! loop track for this meta track
        do k=1,maxmetalen
          if(allmeta(i,k)==-1)exit
          if(trtypes(allmeta(i,k))==10)then
            nkinds(i,1)=nkinds(i,1)+1
          else if(trtypes(allmeta(i,k))==12)then
            nkinds(i,2)=nkinds(i,2)+1
          else if(trtypes(allmeta(i,k))==17)then
            nkinds(i,3)=nkinds(i,3)+1
          else if(trtypes(allmeta(i,k))==18)then
            nkinds(i,4)=nkinds(i,4)+1
          else if(trtypes(allmeta(i,k))==20)then
            nkinds(i,5)=nkinds(i,5)+1
          else if(trtypes(allmeta(i,k))==33)then
            nkinds(i,6)=nkinds(i,6)+1
          else if(trtypes(allmeta(i,k))==34)then
            nkinds(i,7)=nkinds(i,7)+1
          else if(trtypes(allmeta(i,k))==36)then
            nkinds(i,8)=nkinds(i,8)+1
          end if
        end do
      end do

      ! duration
      allocate(metadur(nmeta))
      metadur=-1
      do i=1,nmeta
        mints=ntp+1
        maxts=-1
        do k=1,maxmetalen
          if(allmeta(i,k)==-1)exit
          if(tsclID(alltracks(allmeta(i,k),1))<mints)mints=tsclID(alltracks(allmeta(i,k),1))
          do j=1,maxtrlen
            if(alltracks(allmeta(i,k),j)==-1)exit
            if(tsclID(alltracks(allmeta(i,k),j))>maxts)maxts=tsclID(alltracks(allmeta(i,k),j))
          end do
        end do
        metadur(i)=maxts-mints+1
      end do

      open(unit=1,file="meta_summary.txt",action="write",status="replace")
      write(1,*)" metaID     dur     n10     n12     n17     n18     n20     n33     n34     n36"
      do i=1,nmeta
        write(1,'(10i8)')i,metadur(i),nkinds(i,:)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"==== FINISHED META TRACK STATISTICS ==="
      write(*,*)"======================================="

    end subroutine calcmetatrackstatistics

    subroutine writemetatracks

      use globvar
      use ncdfpars

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=8), allocatable :: dat(:),pdat(:)          ! array for reading float from nc

      write(*,*)"======================================="
      write(*,*)"====== WRITE META TRACKS TO FILE ======"
      write(*,*)"======================================="
      write(*,*)"=== ... to meta.nc ..."
      write(*,*)"---------"

      ! Get initial Information about grid and timesteps of both files
      CALL datainfo(outfile)

      ! Open the dataset 1
      streamID1=streamOpenRead(outfile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if

      ! Set the variable IDs 1
      varID1=0
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)


      !! open new nc file for results
      ! define grid
      gridID2=gridCreate(GRID_GENERIC, nx*ny)
      CALL gridDefXsize(gridID2,nx)
      CALL gridDefYsize(gridID2,ny)
      CALL gridDefXvals(gridID2,xvals)
      CALL gridDefYvals(gridID2,yvals)
      CALL gridDefXunits(gridID2,"m")
      CALL gridDefYunits(gridID2,"m")
      zaxisID2=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID2, level)
      ! define variables
      nmiss2=-999.D0
      vlistID2=vlistCreate()
      varID2=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID2,varID2,"metaID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each meta track")
      CALL vlistDefVarUnits(vlistID2,varID2,"-")
      CALL vlistDefVarMissval(vlistID2,varID2,nmiss2)
      CALL vlistDefVarDatatype(vlistID2,varID2,DATATYPE_INT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite("meta.nc",FILETYPE_NC)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)

      ! allocate data arrays for storage
      allocate(dat(nx*ny),pdat(nx*ny))
      do tsID=0,ntp-1
        ! status to stdout
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1
        end if

        ! Set time step for input and output
        status=streamInqTimestep(streamID1,tsID)
        status=streamDefTimestep(streamID2,tsID)

        ! Read time step from input
        call streamReadVar(streamID1,varID1,dat,nmiss2)

        !cycle if all values are -999
        if(ALL(dat==-999.D0))then
          if(verbose)write(*,*)"NO tracks in timestep:  ",tsID+1
          call streamWriteVar(streamID2,varID2,dat,nmiss2)
          cycle
        end if

        ! now loop dat and replace clIDs with metaIDs
        pdat=-999.D0
        do i=1,nx*ny
          if(dat(i)==-999.D0)then
            cycle
          end if
          if(clmeta(INT(dat(i)))==-1)then
            pdat(i)=-999.D0
            cycle
          end if
          pdat(i)=clmeta(INT(dat(i)))
        end do

        ! write this time step to meta.nc
        call streamWriteVar(streamID2,varID2,pdat,nmiss2)

      end do

      deallocate(dat,pdat)
      CALL gridDestroy(gridID2)
      CALL vlistDestroy(vlistID2)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"===== FINISHED WRITING META TRACKS ===="
      write(*,*)"======================================="

    end subroutine writemetatracks

    subroutine writemetatracksmainstream

      use globvar
      use ncdfpars

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=8), allocatable :: dat(:),pdat(:)          ! array for reading float from nc

      write(*,*)"======================================="
      write(*,*)"====== WRITE MAINSTREAMS TO FILE ======"
      write(*,*)"======================================="
      write(*,*)"=== ... to meta_mainstream.nc ..."
      write(*,*)"---------"

      ! Get initial Information about grid and timesteps of both files
      CALL datainfo(outfile)

      ! Open the dataset 1
      streamID1=streamOpenRead(outfile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if

      ! Set the variable IDs 1
      varID1=0
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)


      !! open new nc file for results
      ! define grid
      gridID2=gridCreate(GRID_GENERIC, nx*ny)
      CALL gridDefXsize(gridID2,nx)
      CALL gridDefYsize(gridID2,ny)
      CALL gridDefXvals(gridID2,xvals)
      CALL gridDefYvals(gridID2,yvals)
      CALL gridDefXunits(gridID2,"m")
      CALL gridDefYunits(gridID2,"m")
      zaxisID2=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID2, level)
      ! define variables
      nmiss2=-999.D0
      vlistID2=vlistCreate()
      varID2=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID2,varID2,"metaID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each meta track")
      CALL vlistDefVarUnits(vlistID2,varID2,"-")
      CALL vlistDefVarMissval(vlistID2,varID2,nmiss2)
      CALL vlistDefVarDatatype(vlistID2,varID2,DATATYPE_INT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite("meta_mainstream.nc",FILETYPE_NC)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)

      ! allocate data arrays for storage
      allocate(dat(nx*ny),pdat(nx*ny))
      do tsID=0,ntp-1
        ! status to stdout
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1
        end if

        ! Set time step for input and output
        status=streamInqTimestep(streamID1,tsID)
        status=streamDefTimestep(streamID2,tsID)

        ! Read time step from input
        call streamReadVar(streamID1,varID1,dat,nmiss2)

        !cycle if all values are -999
        if(ALL(dat==-999.D0))then
          if(verbose)write(*,*)"NO tracks in timestep:  ",tsID+1
          call streamWriteVar(streamID2,varID2,dat,nmiss2)
          cycle
        end if

        ! now loop dat and replace clIDs with metaIDs
        pdat=-999.D0
        do i=1,nx*ny
          if(dat(i)==-999.D0)then
            cycle
          end if
          if(clmetamstr(INT(dat(i)))==-1)then
            pdat(i)=-999.D0
            cycle
          end if
          pdat(i)=clmetamstr(INT(dat(i)))
        end do

        ! write this time step to meta.nc
        call streamWriteVar(streamID2,varID2,pdat,nmiss2)

      end do

      deallocate(dat,pdat)
      CALL gridDestroy(gridID2)
      CALL vlistDestroy(vlistID2)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"===== FINISHED WRITING MAINSTREAMS ===="
      write(*,*)"======================================="

    end subroutine writemetatracksmainstream

end module metatrackstatistics
