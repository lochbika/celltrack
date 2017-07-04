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
      integer                   :: movex,movey             ! the number of gridpoints to move a cell (x and y direction)
      integer, allocatable      :: nCLts(:)

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
      allocate(nCLts(ntp))
      nCLts=0
      maxnIDs=0
      do i=1,globnIDs
        nCLts(tsclID(i))=nCLts(tsclID(i))+1
      end do
      maxnIDs=MAXVAL(nCLts)
      if(MOD(REAL(maxnIDs*4+5),2.0D0)==0.D0)then
        maxnIDs=maxnIDs*4+6
      else
        maxnIDs=maxnIDs*4+5
      end if

      if(globnIDs>25000 .OR. maxnIDs>25000)then
        write(*,*)"=== This may use a lot of RAM! Will allocate: ",maxnIDs*globnIDs*4/1024/1024,"Mb"
      end if

      ! allocate the logical link matrix
      allocate(links(globnIDs,maxnIDs))
      links=.false.
      ! allocate the vector for storing the minclIDs
      allocate(minclIDloc(globnIDs))
      minclIDloc=0

      ! truncate the 2nd dim of the link matrix
      do i=1,globnIDs
        minclIDloc(i)=i-((maxnIDs-1)/2)
        if(minclIDloc(i)<0)minclIDloc(i)=0
      end do
      
      ! now we need to find where in links 2nd dim is i; iclIDloc
      allocate(iclIDloc(globnIDs))
      iclIDloc=1
      do i=1,globnIDs
        do k=minclIDloc(i)+1,minclIDloc(i)+maxnIDs
          if(clIDs(k)==clIDs(i))then
            iclIDloc(i)=k-minclIDloc(i)
            !write(*,*)iclIDloc(i),minclIDloc(i)
            exit
          end if
        end do
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
            allocate(pdat2d(nx,ny),dat2d(nx,ny))
            CALL reshapeT2d(pdat,nx,ny,pdat2d)
            CALL reshapeT2d(dat,nx,ny,dat2d)
            ! allocate temporary storage for advected cells
            allocate(advcell(nx,ny))
            ! move cells in each patch of velocity field
            do x=1,vnx
              do y=1,vny
                ! reset advcell
                advcell=outmissval
                ! calculate the number of grid points to move cells in that patch
                ! movex/movey*(-1) because the shift function does it backwards
                ! do not move if no velocity data is available
                if(uvfield2d(x,y)==outmissval)then 
                  movex=0
                else
                  movex=NINT(uvfield2d(x,y)*tstep/diflon)*(-1)
                end if
                if(vvfield2d(x,y)==outmissval)then
                  movey=0
                else
                  movey=NINT(vvfield2d(x,y)*tstep/diflat)*(-1)
                end if
                ! copy cells from this patch
                do i=1,nx
                  do k=1,ny
                    if(pdat2d(i,k)==outmissval)cycle
                    if(vclxindex(INT(pdat2d(i,k)))==x .AND. vclyindex(INT(pdat2d(i,k)))==y)then
                      advcell(i,k)=pdat2d(i,k)
                    end if
                  end do
                end do
                ! move in x and y direction
                if(periodic)then
                  if(movex.ne.0)advcell=CSHIFT(advcell,SHIFT=movex,DIM=1)
                  if(movey.ne.0)advcell=CSHIFT(advcell,SHIFT=movey,DIM=2)
                else
                  if(movex.ne.0)advcell=EOSHIFT(advcell,SHIFT=movex,BOUNDARY=outmissval,DIM=1)
                  if(movey.ne.0)advcell=EOSHIFT(advcell,SHIFT=movey,BOUNDARY=outmissval,DIM=2)
                end if
                ! do the linking for each patch seperately
                ! now loop all gridpoints
                do i=1,nx
                  do k=1,ny
                    if(advcell(i,k).ne.outmissval .AND. dat2d(i,k).ne.outmissval)then
                      if(verbose)write(*,*)"We have an overlap! cluster ",INT(advcell(i,k))," with ",INT(pdat2d(i,k))
                      j=advcell(i,k)
                      l=dat2d(i,k)
                      ! backward linking
                      links(l,j-minclIDloc(l))=.true.
                      ! forward linking
                      links(j,l-minclIDloc(j))=.true.
                    end if
                  end do
                end do
              end do
            end do
            deallocate(pdat2d,dat2d,advcell)
          else
            ! now loop all gridpoints without advection correction
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
          end if
          deallocate(pdat)
        end if
        deallocate(dat)
        if(advcor .AND. adviter>1)then
          deallocate(uvfield2d,vvfield2d)
        end if
      end do

      CALL streamClose(streamID2)
      if(advcor .AND. adviter>1)CALL streamClose(streamID3)

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
