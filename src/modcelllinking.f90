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
      real(kind=8), allocatable :: dat(:),pdat(:)          ! array for reading float from nc
      real(kind=8), allocatable :: pdat2d(:,:)  ! array for doing the clustering

      write(*,*)"======================================="
      write(*,*)"=========== STARTED LINKING ==========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Searching for fw/bw links betw cells ..."
      write(*,*)"---------"
    
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
        if(tsID==1)then
          minclIDloc(i)=0
          cycle
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
    
      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
          write(*,*)"Processing timestep: ",tsID+1
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
    
        if(tsID.ne.0)then
          ! load previous timestep
          allocate(pdat(nx*ny))
          ! Set time step for input file
          status=streamInqTimestep(streamID1,tsID-1)
          ! Read time step from input
          call streamReadVar(streamID1,varID1,pdat,nmiss2)
          ! now loop all gridpoints
          do i=1,nx*ny
            if(dat(i).ne.-999.D0 .AND. pdat(i).ne.-999.D0)then
              if(verbose)write(*,*)"BACKWARD::We have an overlap! cluster ",dat(i)," with ",pdat(i)
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
      end do
    
      CALL streamClose(streamID1)
    
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
    
      if(lout)then
    
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
