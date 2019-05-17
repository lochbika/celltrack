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
module cellshape

  use globvar

  implicit none

  contains
    subroutine calccellshape

      use ncdfpars
      use globvar

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=stdfloattype), allocatable :: dat(:),dat2d(:,:) ! array for reading cellIDs from nc
      real(kind=stdfloattype), allocatable :: coords(:,:,:)     ! array for holding all cells coordinates
      integer, allocatable :: cellcnt(:)             ! array for holding a counter for each cell
      real(kind=stdfloattype), allocatable  :: axisLen(:,:)
      real(kind=stdfloattype), allocatable  :: rot(:)           ! rotation angles in radians

      real(kind=stdfloattype) :: distxo,distyo ! the direct distance if cells do not cross boundaries
      real(kind=stdfloattype) :: distxp,distyp ! the distance between if cells touch boundaries

      integer :: minclID,maxclID,maxarea
      real(kind=stdfloattype) :: maxLen
      real(kind=stdfloattype), allocatable  :: tmpcoords(:,:)

      write(*,*)"======================================="
      write(*,*)"============= CELL SHAPE =============="
      write(*,*)"======================================="

      ! the rotation angles
      allocate(rot(361))
      rot(1)=0
      do i=2,361
        rot(i)=(2*pi)*(i-1)/360
        !write(*,*)rot(i)
      end do

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

      ! allocate array for the length of the minor and major axis
      allocate(axisLen(globnIDs,3))
      axisLen=0

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
        allocate(cellcnt(minclID:maxclID))
        coords=-1
        cellcnt=0

        ! reshape to 2D
        allocate(dat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        deallocate(dat)
        
        ! assign coordinates
        do x=1,nx
          do y=1,ny
            if(dat2d(x,y).ne.outmissval)then
              cellcnt(INT(dat2d(x,y)))=cellcnt(INT(dat2d(x,y)))+1
              coords(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),1)=x
              coords(INT(dat2d(x,y)),cellcnt(INT(dat2d(x,y))),2)=y
            end if
          end do
        end do
        deallocate(dat2d)
        
        ! center each cell on the center of mass
        ! if cell touches domain boundaries and we have periodic boundary conditions
        ! -> use the shorter distance like in modadvectioncorrection
        do clID=minclID,maxclID
          do i=1,cellcnt(clID)
            if(periodic)then
                ! ok, now it's getting tricky because we don't know whether the cell
                ! wraps around the boundaries.
                ! At this step I calculate two kinds of distances
                ! 1. the direct distance
                !    this assumes the center of mass and the grid point are not across boundaries 
                ! 2. the distance between center of mass and the projection of the grid point
                !    this assumes that both points are closer across boundaries
                distxo=coords(clID,i,1)-clcmass(clID,1) ! 1 in x dir
                distyo=coords(clID,i,2)-clcmass(clID,2) ! 1 in y dir
                ! 2 in x direction
                if(clcmass(clID,1)>=coords(clID,i,1))then
                  distxp=coords(clID,i,1)+nx-clcmass(clID,1)
                else
                  distxp=(clcmass(clID,1)+nx-coords(clID,i,1))*(-1.0D0)
                end if
                ! 2 in y direction
                if(clcmass(clID,2)>=coords(clID,i,2))then
                  distyp=coords(clID,i,2)+ny-clcmass(clID,2)
                else
                  distyp=(clcmass(clID,2)+ny-coords(clID,i,2))*(-1.0D0)
                end if
                if(verbose)write(*,*)"distances are x,y"
                if(verbose)write(*,*)"  ",distxo,distyo
                if(verbose)write(*,*)"projected distances are x,y"
                if(verbose)write(*,*)"  ",distxp,distyp
                if( abs(distxo) .gt. abs(distxp) )then
                  ! this means both points are closer across boundaries
                  coords(clID,i,1)=distxp
                  if(verbose)write(*,*)"Cell ",clID," crosses the x boundary"
                else
                  ! no boundary wrapping: just do the normal calculation
                  coords(clID,i,1)=distxo
                end if
                if( abs(distyo) .gt. abs(distyp) )then
                  ! this means both points are closer across boundaries
                  coords(clID,i,2)=distyp
                  if(verbose)write(*,*)"Cell ",clID," crosses the y boundary"
                else
                  ! no boundary wrapping: just do the normal calculation
                  coords(clID,i,2)=distyo
                end if
            else
              coords(clID,i,1)=coords(clID,i,1)-clcmass(clID,1)
              coords(clID,i,2)=coords(clID,i,2)-clcmass(clID,2)
            end if
          end do
        end do

        ! rotate each cell in steps of 1 degree
        ! x' = x*cos(theta) + y*sin(theta)
        ! y' = -1 * x*sin(theta) + y*cos(theta)
        do clID=minclID,maxclID
          ! rotate and save coordinates of the current rotation
          maxLen=0
          do j=1,361
            allocate(tmpcoords(cellcnt(clID),2))
            do i=1,cellcnt(clID)
              tmpcoords(i,1)=( coords(clID,i,1) * COS(rot(j)) ) + ( coords(clID,i,2) * SIN(rot(j)) )
              tmpcoords(i,2)=( -1 * coords(clID,i,1) * SIN(rot(j)) ) + ( coords(clID,i,2) * COS(rot(j)) )
              !if(j==180)write(*,*)clID,tmpcoords(i,:),coords(clID,i,:),rot(j)
            end do
            ! for this rotation find the longest distance in x direction and the corresponding length in y direction
            if(((MAXVAL(tmpcoords(:,1))+.5) - (MINVAL(tmpcoords(:,1))-.5))>maxLen)then
              maxLen=(MAXVAL(tmpcoords(:,1))+.5) - (MINVAL(tmpcoords(:,1))-.5)
              axisLen(clID,1)=MAX(maxlen,1.D0)
              axisLen(clID,2)=MAX((MAXVAL(tmpcoords(:,2))+.5) - (MINVAL(tmpcoords(:,2))-.5),1.D0)
              axisLen(clID,3)=j-1
            end if
            deallocate(tmpcoords)
          end do
        end do

        deallocate(coords,cellcnt)
      end do

      CALL streamClose(streamID2)

      ! now write it to a file
      open(unit=1,file="cell_shape.txt",action="write",status="replace")
      write(1,*)"     clID           majorLen           minorLen               angle"
      do i=1,globnIDs
        write(1,'(1i10,3f19.6)')clIDs(i),axisLen(i,1),axisLen(i,2),axisLen(i,3)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"========= FINISHED CELL SHAPE ========="
      write(*,*)"======================================="

    end subroutine calccellshape

end module cellshape
