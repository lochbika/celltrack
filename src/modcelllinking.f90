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
module celllinking

  use globvar

  implicit none

  contains
    subroutine linking

      use globvar
      use ncdfpars

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=8), allocatable :: dat(:),pdat(:)          ! arrays for reading float from nc
      real(kind=8), allocatable :: dat2d(:,:),pdat2d(:,:)  ! arrays for the advection correction
      real(kind=8), allocatable :: advcell(:,:)            ! temporary array for advected cells
      real(kind=8), allocatable :: preadvcell(:,:)         ! temporary array for cells of patches
      integer :: movex,movey                               ! the number of gridpoints to move a cell (x and y direction)

      write(*,*)"======================================="
      write(*,*)"=========== STARTED LINKING ==========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Searching for fw/bw links betw cells ..."
      write(*,*)"---------"

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

      ! find th maximum number of cells per timestep
      ! to reduce 2nd dimension of the logical link matrix
      maxnIDs=0
      tp=0
      k=0
      do i=1,globnIDs
        j=tsclID(i)
        if(j==k)tp=tp+1
        if(j>k)then
          if(maxnIDs<tp)maxnIDs=tp
          k=j
          tp=0
        end if
      end do
      maxnIDs=maxnIDs*3+5

      if(globnIDs>25000 .OR. maxnIDs>25000)then
        write(*,*)"=== This may use a lot of RAM! Will allocate: ",maxnIDs*globnIDs*4/1024/1024,"Mb"
      end if

      ! allocate the logical link matrix
      allocate(links(globnIDs,maxnIDs))
      links=.false.
      ! allocate the vector for storing the minclIDs
      allocate(minclIDloc(globnIDs))
      minclIDloc=-1

      ! truncate the 2nd dim of the link matrix
      do i=1,globnIDs
        tsID=tsclID(i)
        if(tsID==1 .OR. i<(maxnIDs-5)/3)then
          minclIDloc(i)=0
        else
          do k=1,i
            if(tsclID(k)==tsID-1 .AND. tsclID(k)>1)then
              minclIDloc(i)=k-1
              exit
            end if
            ! backup: if there are no cells in the previous timestep
            if(tsclID(k)==tsID .AND. tsclID(k)>1)then
              minclIDloc(i)=k-1
              exit
            end if
            ! backup: if there are no cells in the current timestep
            if(tsclID(k)==tsID+1 .AND. tsclID(k)>1)then
              minclIDloc(i)=-1
              exit
            end if
          end do
        end if
      end do

      ! if we do advection correction... read the vfile now
      if(advcor .AND. adviter>1)then
        write(vfile,'(A7,I0.3,A3)')"vfield_",adviter-1,".nc"

        write(*,*)"======================================="
        write(*,*)"=== Opening file with velocity fields ",TRIM(vfile)
        write(*,*)"---------"

        streamID3=streamOpenRead(TRIM(vfile))
        if(streamID3<0)then
           write(*,*)cdiStringError(streamID3)
           stop
        end if
        vuID=0
        vvID=1
        vlistID3=streamInqVlist(streamID3)
        gridID3=vlistInqVarGrid(vlistID2,vuID)
        taxisID3=vlistInqTaxis(vlistID3)
        zaxisID3=vlistInqVarZaxis(vlistID3,vuID)
        outmissval=vlistInqVarMissval(vlistID3,vuID)
      end if

      ! do the linking per time step
      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1
        end if

        ! Allocate arrays for data storage
        allocate(dat(nx*ny))
        if(advcor .AND. adviter>1)then
          allocate(uvfield(vnx*vny),vvfield(vnx*vny))
          allocate(uvfield2d(vnx,vny),vvfield2d(vnx,vny))
        end if

        ! Set time step for input files
        status=streamInqTimestep(streamID2,tsID)
        if(advcor .AND. adviter>1)status=streamInqTimestep(streamID3,tsID)

        ! Read time step from input
        call streamReadVar(streamID2,varID2,dat,nmiss2)
        if(advcor .AND. adviter>1)then
          CALL streamReadVar(streamID3,vuID,uvfield,nmiss3)
          CALL streamReadVar(streamID3,vvID,vvfield,nmiss3)
          CALL reshapeT2d(uvfield,vnx,vny,uvfield2d)
          CALL reshapeT2d(vvfield,vnx,vny,vvfield2d)
          deallocate(uvfield,vvfield)
        end if

        !cycle if all values are -999
        if(nmiss2==nx*ny)then
          if(verbose)write(*,*)"NO clusters in timestep:  ",tsID+1
          deallocate(dat)
          if(advcor .AND. adviter>1)then
            deallocate(uvfield2d,vvfield2d)
          end if
          cycle
        end if

        if(tsID.ne.0)then
          ! load previous timestep
          allocate(pdat(nx*ny))
          ! Set time step for input file
          status=streamInqTimestep(streamID2,tsID-1)
          ! Read time step from input
          call streamReadVar(streamID2,varID2,pdat,nmiss2)

          ! ADVECTION CORRECTION: now advect previous timestep with velocity field
          if(advcor .AND. adviter>1)then
            ! reshape to 2d
            allocate(pdat2d(nx,ny))
            CALL reshapeT2d(pdat,nx,ny,pdat2d)
            ! allocate temporary storage for advected cells
            allocate(advcell(nx,ny))
            advcell=outmissval
            allocate(preadvcell(nx,ny))
            ! move cells in each patch of velocity field
            do x=1,vnx
              do y=1,vny
                ! cycle if no velocity data is available
                if(uvfield2d(x,y)==outmissval .AND. vvfield2d(x,y)==outmissval)cycle
                ! calculate the number of grid points to move cells in that patch
                ! movex/movey*(-1) because the shift function does it backwards
                movex=NINT(uvfield2d(x,y)*tstep/diflon)*(-1)
                movey=NINT(vvfield2d(x,y)*tstep/diflat)*(-1)
                ! now reset preadvcell and copy cells inside this patch to it
                preadvcell=outmissval
                do i=1,nx
                  do k=1,ny
                    if(pdat2d(i,k)==outmissval)cycle
                    if(vclxindex(INT(pdat2d(i,k)))==x .AND. vclyindex(INT(pdat2d(i,k)))==y)then
                      preadvcell(i,k)=pdat2d(i,k)
                    end if
                  end do
                end do
                ! move in x and y direction
                if(movex.ne.0)then
                  preadvcell=EOSHIFT(preadvcell,SHIFT=movex,BOUNDARY=outmissval,DIM=1)
                end if
                if(movey.ne.0)then
                  preadvcell=EOSHIFT(preadvcell,SHIFT=movey,BOUNDARY=outmissval,DIM=2)
                end if
                ! now copy the advected cells from this patch to advcell
                do i=1,nx
                  do k=1,ny
                    if(preadvcell(i,k)==outmissval)cycle
                    if(vclxindex(INT(preadvcell(i,k)))==x .AND. vclyindex(INT(preadvcell(i,k)))==y)then
                      advcell(i,k)=preadvcell(i,k)
                    end if
                  end do
                end do
              end do
            end do
            ! now bring it back home
            CALL reshapeF2d(advcell,nx,ny,pdat)
            deallocate(pdat2d,preadvcell,advcell)
          end if

          ! now loop all gridpoints
          do i=1,nx*ny
            if(dat(i).ne.outmissval .AND. pdat(i).ne.outmissval)then
              if(verbose)write(*,*)"We have an overlap! cluster ",INT(dat(i))," with ",INT(pdat(i))
              k=dat(i)
              j=pdat(i)
              ! backward linking
              links(k,j-minclIDloc(k))=.true.
              ! forward linking
              links(j,k-minclIDloc(j))=.true.
            end if
          end do
          deallocate(pdat)
        end if
        deallocate(dat)
        if(advcor .AND. adviter>1)then
          deallocate(uvfield2d,vvfield2d)
        end if
      end do

      CALL streamClose(streamID2)
      if(advcor .AND. adviter>1)CALL streamClose(streamID3)

      ! now we need to find where in links 2nd dim is i; iclIDloc
      allocate(iclIDloc(globnIDs))
      iclIDloc=-1
      do i=1,globnIDs
        if(tsclID(i).ne.0 .AND. minclIDloc(i).ne.-1)then
          do k=minclIDloc(i)+1,minclIDloc(i)+maxnIDs+1
            if(k>globnIDs)exit
            if(clIDs(k)==clIDs(i))then
              iclIDloc(i)=k-minclIDloc(i)
              exit
            end if
          end do
        end if
      end do

      if(lout .AND. adviter==nadviter+1)then

        write(*,*)"======================================="
        write(*,*)"=== writing links to file cell_links.txt ..."
        write(*,*)"---------"

        ! write the link matrix to file
        open(unit=1,file="cell_links.txt",action="write",status="replace")
        write(1,*)"stIDloc   stID      .......        "
        do i=1,globnIDs
          write(1,*)minclIDloc(i),iclIDloc(i),clIDs(i),links(i,:)
        end do
        close(unit=1)

      end if

      write(*,*)"======================================="
      write(*,*)"========== FINISHED LINKING ==========="
      write(*,*)"======================================="

  end subroutine linking

end module celllinking
