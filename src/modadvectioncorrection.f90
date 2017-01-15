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
      real(kind=8), allocatable :: smpsize2d(:,:),smpsize(:)   ! sample size for each gridpoint on the velocity field
      real(kind=8) :: mindist,cdist

      write(*,*)"======================================="
      write(*,*)"===== START ADVECTION CORRECTION ======"
      write(*,*)"======================================="

      allocate(vclxindex(globnIDs),vclyindex(globnIDs))
      allocate(vclx(globnIDs),vcly(globnIDs))

      ! we can use the gathered information to coarse grain the grid and open a new dataset
      vnx=nx/coarsex
      vny=ny/coarsey
      allocate(vxvals(vnx),vyvals(vny))
      do x=1,vnx
        vxvals(x)=(xvals(0)-diflon/2) + (diflon*coarsex*x) - (diflon*coarsex/2)
      end do
      do y=1,vny
        vyvals(y)=(yvals(0)-diflat/2) + (diflat*coarsey*y) - (diflat*coarsey/2)
      end do
      
      write(*,*)"======================================="
      write(*,*)"=== GRID FOR ADVECTION CORRECTION:"
      write(*,'(A,1i12)')" NX      : ",vnx
      write(*,'(A,1f12.2)')" MIN X   : ",vxvals(1)
      write(*,'(A,1f12.2)')" MAX X   : ",vxvals(vnx)
      write(*,'(A,1f12.2)')" DIF X   : ",(vxvals(1)-xvals(0))*2
      write(*,'(A,1a12)')" Unit    : ",trim(xunit)
      write(*,'(A,1i12)')" NY      : ",vny
      write(*,'(A,1f12.2)')" MIN Y   : ",vyvals(1)
      write(*,'(A,1f12.2)')" MAX Y   : ",vyvals(vny)
      write(*,'(A,1f12.2)')" DIF Y   : ",(vyvals(1)-yvals(0))*2
      write(*,'(A,1a12)')" Unit    : ",trim(yunit)
      write(*,*)"---------"

      ! find the nearest gridpoint on the velocity grid for all cells
      do clID=1,globnIDs
        mindist=HUGE(mindist)
        do x=1,vnx
          cdist=abs(vxvals(x)-(wclcmass(clID,1)*diflon+xvals(0)))
          if(cdist<mindist)then
            vclxindex(clID)=x
            mindist=cdist
          end if
        end do
        mindist=HUGE(mindist)
        do y=1,vny
          cdist=abs(vyvals(y)-(wclcmass(clID,2)*diflat+yvals(0)))
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
        varID1=ivar
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
        CALL vlistDefVarDatatype(vlistID2,vuID,DATATYPE_FLT64)
        vvID=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
        CALL vlistDefVarName(vlistID2,vvID,"v")
        CALL vlistDefVarLongname(vlistID2,vvID,"derived wind speed in y direction")
        CALL vlistDefVarUnits(vlistID2,vvID,"m/s")
        CALL vlistDefVarMissval(vlistID2,vvID,outmissval)
        CALL vlistDefVarDatatype(vlistID2,vvID,DATATYPE_FLT64)
        ssizeID=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
        CALL vlistDefVarName(vlistID2,ssizeID,"sample_size")
        CALL vlistDefVarLongname(vlistID2,ssizeID,"number of cells in grid area")
        CALL vlistDefVarUnits(vlistID2,ssizeID,"-")
        CALL vlistDefVarMissval(vlistID2,ssizeID,outmissval)
        CALL vlistDefVarDatatype(vlistID2,ssizeID,DATATYPE_INT32)
        ! copy time axis from input
        taxisID2=vlistInqTaxis(vlistID1)
        call vlistDefTaxis(vlistID2,taxisID2)
       
        ! Open the dataset for writing
        write(vfile,'(A7,I0.3,A3)')"vfield_",adviter,".nc"
        streamID2=streamOpenWrite(TRIM(vfile),FILETYPE_NC)
        if(streamID2<0)then
           write(*,*)cdiStringError(streamID2)
           stop
        end if
        
        write(*,*)"======================================="
        write(*,*)"=== Calc velocity field and write to ",TRIM(vfile),"..."
        write(*,*)"---------"

        ! Assign variables to dataset
        call streamDefVList(streamID2,vlistID2)

        ! now calculate each cells velocity
        vclx=outmissval
        vcly=outmissval
        do clID=1,globnIDs
          if(tsclID(clID).ne.1 .AND. nbw(clID)==1)then
            ! find the cell which is connected backwards
            do i=1,iclIDloc(clID)
              if(links(clID,i))then
                selCL=i+minclIDloc(clID)
                exit
              end if
            end do
            if(nfw(selCL)==1)then
              vclx(clID)=(wclcmass(clID,1)-wclcmass(selCL,1))*diflon/tstep
              vcly(clID)=(wclcmass(clID,2)-wclcmass(selCL,2))*diflat/tstep
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
        deallocate(links,minclIDloc,iclIDloc,nbw,nfw)

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
