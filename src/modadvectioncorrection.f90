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
module advectioncorrection

  use globvar

  implicit none

  contains
    subroutine doadvectioncorrection

      ! modules
      use globvar
      use ncdfpars
      use celllinking
      use linkstatistics

      implicit none

      include 'cdi.inc'

      ! variables and arrays
      integer :: selCL
      real(kind=stdfloattype), allocatable :: smpsize2d(:,:),smpsize(:)   ! sample size for each gridpoint on the velocity field
      real(kind=stdfloattype) :: mindist,cdist,vxres,vyres
      real(kind=stdfloattype) :: distxo,distyo ! the direct distance if cells do not cross boundaries
      real(kind=stdfloattype) :: distxp,distyp ! the distance between if cells crossed boundaries

      write(*,*)"======================================="
      write(*,*)"===== START ADVECTION CORRECTION ======"
      write(*,*)"======================================="

      allocate(vclxindex(globnIDs),vclyindex(globnIDs))
      allocate(vclx(globnIDs),vcly(globnIDs))

      ! we can use the gathered information to coarse grain the grid and open a new dataset
      vnx=nx/coarsex
      vny=ny/coarsey
      vxres=( (xvals(nx)+diflon/2) - (xvals(1)-diflon/2) ) / vnx
      vyres=( (yvals(ny)+diflat/2) - (yvals(1)-diflat/2) ) / vny
      allocate(vxvals(vnx),vyvals(vny))
      vxvals(1)= (xvals(1)-diflon/2) + (vxres/2)
      vyvals(1)= (yvals(1)-diflat/2) + (vyres/2)
      if(vnx>1)then
        do x=2,vnx
          vxvals(x)=vxvals(1) + (x-1)*vxres
        end do
      end if
      if(vny>1)then
        do y=2,vny
          vyvals(y)=vyvals(1) + (y-1)*vyres
        end do
      end if

      write(*,*)"======================================="
      write(*,*)"=== GRID FOR ADVECTION CORRECTION:"
      write(*,'(A,1i12)')" NX      : ",vnx
      write(*,'(A,1f12.2)')" MIN X   : ",vxvals(1)
      write(*,'(A,1f12.2)')" MAX X   : ",vxvals(vnx)
      write(*,'(A,1f12.2)')" DIF X   : ",vxres
      write(*,'(A,1a12)')" Unit    : ",trim(xunit)
      write(*,'(A,1i12)')" NY      : ",vny
      write(*,'(A,1f12.2)')" MIN Y   : ",vyvals(1)
      write(*,'(A,1f12.2)')" MAX Y   : ",vyvals(vny)
      write(*,'(A,1f12.2)')" DIF Y   : ",vyres
      write(*,'(A,1a12)')" Unit    : ",trim(yunit)
      write(*,*)"---------"

      ! find the nearest gridpoint on the velocity grid for all cells
      vclxindex=-1
      vclyindex=-1
      do clID=1,globnIDs
        mindist=HUGE(mindist)
        do x=1,vnx
          cdist=abs(vxvals(x)-wclcmass(clID,1))
          if(cdist<mindist)then
            vclxindex(clID)=x
            mindist=cdist
          end if
        end do
        mindist=HUGE(mindist)
        do y=1,vny
          cdist=abs(vyvals(y)-wclcmass(clID,2))
          if(cdist<mindist)then
            vclyindex(clID)=y
            mindist=cdist
          end if
        end do
      end do

      do adviter=1,nadviter

        write(*,*)"======================================="
        write(*,*)"=== This is iteration ",adviter
        write(*,*)"---------"

        ! now we do the linking and the statistics
        CALL linking()
        CALL calclinkstatistics()

        ! Open the dataset 1
        streamID1=streamOpenRead(ifile)
        if(streamID1<0)then
           write(*,*)cdiStringError(streamID1)
           stop
        end if

        ! Set the variable IDs 1
        varID1=getVarIDbyName(ifile,ivar)
        vlistID1=streamInqVlist(streamID1)
        gridID1=vlistInqVarGrid(vlistID1,varID1)
        taxisID1=vlistInqTaxis(vlistID1)
        zaxisID1=vlistInqVarZaxis(vlistID1,varID1)

        !! open new nc file for results
        ! define grid
        gridID2=gridCreate(GRID_GENERIC, vnx*vny)
        CALL gridDefXsize(gridID2,vnx)
        CALL gridDefYsize(gridID2,vny)
        CALL gridDefXvals(gridID2,vxvals)
        CALL gridDefYvals(gridID2,vyvals)
        CALL gridDefXunits(gridID2,TRIM(xunit))
        CALL gridDefYunits(gridID2,TRIM(yunit))
        zaxisID2=zaxisCreate(ZAXIS_GENERIC, 1)
        CALL zaxisDefLevels(zaxisID2, level)
        ! define variables
        vlistID2=vlistCreate()
        vuID=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
        CALL vlistDefVarName(vlistID2,vuID,"u")
        CALL vlistDefVarLongname(vlistID2,vuID,"derived wind speed in x direction")
        CALL vlistDefVarUnits(vlistID2,vuID,"m/s")
        CALL vlistDefVarMissval(vlistID2,vuID,outmissval)
        CALL vlistDefVarDatatype(vlistID2,vuID,CDI_DATATYPE_FLT64)
        vvID=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
        CALL vlistDefVarName(vlistID2,vvID,"v")
        CALL vlistDefVarLongname(vlistID2,vvID,"derived wind speed in y direction")
        CALL vlistDefVarUnits(vlistID2,vvID,"m/s")
        CALL vlistDefVarMissval(vlistID2,vvID,outmissval)
        CALL vlistDefVarDatatype(vlistID2,vvID,CDI_DATATYPE_FLT64)
        ssizeID=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
        CALL vlistDefVarName(vlistID2,ssizeID,"sample_size")
        CALL vlistDefVarLongname(vlistID2,ssizeID,"number of cells in grid area")
        CALL vlistDefVarUnits(vlistID2,ssizeID,"-")
        CALL vlistDefVarMissval(vlistID2,ssizeID,outmissval)
        CALL vlistDefVarDatatype(vlistID2,ssizeID,CDI_DATATYPE_INT32)
        ! copy time axis from input
        taxisID2=vlistInqTaxis(vlistID1)
        call vlistDefTaxis(vlistID2,taxisID2)

        ! Open the dataset for writing
        write(vfile,'(A7,I0.3,A3)')"vfield_",adviter,".nc"
        streamID2=streamOpenWrite(TRIM(vfile),CDI_FILETYPE_NC4)
        if(streamID2<0)then
           write(*,*)cdiStringError(streamID2)
           stop
        end if

        write(*,*)"======================================="
        write(*,*)"=== Calc velocity field and write to ",TRIM(vfile),"..."
        write(*,*)"---------"

        ! set netCDF4 compression
        CALL streamDefCompType(streamID1,CDI_COMPRESS_ZIP)
        CALL streamDefCompLevel(streamID1, 6)
        ! Assign variables to dataset
        call streamDefVList(streamID2,vlistID2)

        ! now calculate each cells velocity
        vclx=outmissval
        vcly=outmissval
        do clID=1,globnIDs
          if(touchb(clID))cycle
          if(tsclID(clID).ne.1 .AND. nbw(clID)==1)then
            ! find the cell which is connected backwards
            do i=1,nlinks(clID)
              if(ltype(clID,i)==1)then
                selCL=links(clID,i)
                exit
              end if
            end do
            if(nfw(selCL)==1)then
              if(verbose)write(*,*)"Calculating distance and velocity between cell ",clID," and ",selCL
              if(periodic)then
                ! ok, now it's getting tricky because we don't know whether selCL
                ! crossed the boundaries to become clID.
                ! At this step I calculate two kinds of distances between the two cells
                ! 1. the direct distance: wclcmass(clID,1)-wclcmass(selCL,1)
                !    this assumes selCL did not cross the boundaries
                ! 2. the distance between clID and the projection of selCL
                !    this assumes that selCL crossed the boundaries
                distxo=wclcmassgrd(clID,1)-wclcmassgrd(selCL,1) ! 1 in x direction
                distyo=wclcmassgrd(clID,2)-wclcmassgrd(selCL,2) ! 1 in y direction
                ! 2 in x direction
                if(wclcmassgrd(clID,1)>=wclcmassgrd(selCL,1))then
                  distxp=wclcmassgrd(clID,1)-wclcmassgrd(selCL,1)+nx
                else
                  distxp=wclcmassgrd(clID,1)+nx-wclcmassgrd(selCL,1)
                end if
                ! 2 in y direction
                if(wclcmassgrd(clID,2)>=wclcmassgrd(selCL,2))then
                  distyp=wclcmassgrd(clID,2)-wclcmassgrd(selCL,2)+ny
                else
                  distyp=wclcmassgrd(clID,2)+ny-wclcmassgrd(selCL,2)
                end if
                if(verbose)write(*,*)"distances are x,y"
                if(verbose)write(*,*)"  ",distxo,distyo
                if(verbose)write(*,*)"projected distances are x,y"
                if(verbose)write(*,*)"  ",distxp,distyp
                if( abs(distxo) .gt. abs(distxp) )then
                  ! this means the cell crossed the boundaries
                  vclx(clID)=distxp*diflon/tstep
                  if(verbose)write(*,*)"Cell ",clID," crosses the x boundary"
                else
                  ! no boundary crossing: just do the normal calculation
                  vclx(clID)=distxo*diflon/tstep
                end if
                if( abs(distyo) .gt. abs(distyp) )then
                  ! this means the cell crossed the boundaries
                  vcly(clID)=distyp*diflat/tstep
                  if(verbose)write(*,*)"Cell ",clID," crosses the y boundary"
                else
                  ! no boundary crossing: just do the normal calculation
                  vcly(clID)=distyo*diflat/tstep
                end if
              else
                vclx(clID)=(wclcmass(clID,1)-wclcmass(selCL,1))/tstep
                vcly(clID)=(wclcmass(clID,2)-wclcmass(selCL,2))/tstep
              end if
            end if
          end if
        end do

        ! calculate average velocities on the grid for each time step
        do tsID=0,(ntp-1)
          allocate(uvfield2d(vnx,vny),vvfield2d(vnx,vny),smpsize2d(vnx,vny))
          uvfield2d=0
          vvfield2d=0
          smpsize2d=0

          do clID=1,globnIDs
            if(tsclID(clID)>tsID+1)exit
            if(tsclID(clID)==tsID+1 .AND. vclx(clID).ne.outmissval .AND. vcly(clID).ne.outmissval)then
              if(sqrt(vclx(clID)**2 + vcly(clID)**2)>maxvel)then
                if(verbose)write(*,*)"Skipping cell ",clID," because it has a really high velocity: ", &
                 & sqrt(vclx(clID)**2 + vcly(clID)**2)
                cycle ! cycle if this cell has a unrealistically high velocity
              end if
              uvfield2d(vclxindex(clID),vclyindex(clID)) = uvfield2d(vclxindex(clID),vclyindex(clID)) + vclx(clID)
              vvfield2d(vclxindex(clID),vclyindex(clID)) = vvfield2d(vclxindex(clID),vclyindex(clID)) + vcly(clID)
              smpsize2d(vclxindex(clID),vclyindex(clID)) = smpsize2d(vclxindex(clID),vclyindex(clID)) + 1
            end if
          end do

          ! average and set 0 sampled gridpoints to missing value
          WHERE(smpsize2d.ne.0)uvfield2d=uvfield2d/smpsize2d
          WHERE(smpsize2d.ne.0)vvfield2d=vvfield2d/smpsize2d
          WHERE(smpsize2d==0)uvfield2d=outmissval
          WHERE(smpsize2d==0)vvfield2d=outmissval

          ! reshape to 2D
          allocate(uvfield(vnx*vny),vvfield(vnx*vny),smpsize(vnx*vny))
          CALL reshapeF2d(uvfield2d,vnx,vny,uvfield)
          CALL reshapeF2d(vvfield2d,vnx,vny,vvfield)
          CALL reshapeF2d(smpsize2d,vnx,vny,smpsize)
          deallocate(uvfield2d,vvfield2d,smpsize2d)

          ! now write to vfile
          status=streamDefTimestep(streamID2,tsID)
          CALL streamWriteVar(streamID2,vuID,uvfield,nmiss2)
          CALL streamWriteVar(streamID2,vvID,vvfield,nmiss2)
          CALL streamWriteVar(streamID2,ssizeID,smpsize,nmiss2)

          deallocate(uvfield,vvfield,smpsize)

        end do

        ! deallocate some arrays to rerun the linking part
        deallocate(links,nlinks,clink,ltype,nbw,nfw)

        ! close input and output
        CALL gridDestroy(gridID2)
        CALL vlistDestroy(vlistID2)
        CALL streamClose(streamID2)
        CALL streamClose(streamID1)

      end do

      ! set adviter that the latest vfield file will be used later
      adviter=nadviter+1

      write(*,*)"======================================="
      write(*,*)"==== FINISHED ADVECTION CORRECTION ===="
      write(*,*)"======================================="

    end subroutine doadvectioncorrection

end module advectioncorrection
