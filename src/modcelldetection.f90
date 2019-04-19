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
module celldetection
  
  use globvar

  implicit none
  
  
  contains 
    subroutine docelldetection
      
      use globvar
      use ncdfpars
      
      implicit none
      
      include 'cdi.inc'      
      
      integer :: nIDs,globID
      real(kind=stdfloattype), allocatable :: cl(:,:)
      
      ! data arrays
      real(kind=stdfloattype), allocatable :: dat(:)          ! array for reading float from nc
      real(kind=stdfloattype), allocatable :: dat2d(:,:)      ! array for doing the clustering
      
      globnIDs=0
      globID=1
      
      write(*,*)"======================================="
      write(*,*)"======== START CELL DETECTION ========="
      write(*,*)"======================================="
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to input file..."
    
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
      
      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(outfile)
      write(*,*)"---------"
    
      !! open new nc file for results
      ! define grid
      gridID2=gridDuplicate(gridID1)
      !CALL gridDefXsize(gridID2,nx)
      !CALL gridDefYsize(gridID2,ny)
      !CALL gridDefXvals(gridID2,xvals)
      !CALL gridDefYvals(gridID2,yvals)
      !CALL gridDefXunits(gridID2,TRIM(xunit))
      !CALL gridDefYunits(gridID2,TRIM(yunit))
      zaxisID2=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID2, level)
      ! define variables
      vlistID2=vlistCreate()
      varID2=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID2,varID2,"cellID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each cell")
      CALL vlistDefVarUnits(vlistID2,varID2,"-")
      CALL vlistDefVarMissval(vlistID2,varID2,outmissval)
      CALL vlistDefVarDatatype(vlistID2,varID2,CDI_DATATYPE_INT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite(trim(outfile),CDI_FILETYPE_NC4)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! set netCDF4 compression
      CALL streamDefCompType(streamID2,CDI_COMPRESS_ZIP)
      CALL streamDefCompLevel(streamID2, 6)
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)
      
      ! save which time steps are completely filled with NA
      allocate(tsALLna(ntp))
      tsALLna=.false.

      ! Get data
      write(*,*)"======================================="
      write(*,*)"=== Find continous cells ..."
      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1 .OR. verbose)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if
    
        ! Allocate arrays for data storage
        allocate(dat(nx*ny))
    
        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)

        ! set time step for output file
        status=streamDefTimestep(streamID2,tsID)

        ! Read time step from input
        call streamReadVarSlice(streamID1,varID1,levelID,dat,nmiss1)
        
        ! cycle if field contains only missing values; but write it to output
        if(nmiss1==nx*ny)then
          nmiss2=nmiss1
          dat=outmissval
          tsALLna(tsID+1)=.true.
          CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
          deallocate(dat)
          cycle
        end if
        
        ! reshape array
        allocate(dat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        deallocate(dat)
    
        ! cluster the frame
        ! it is very important that at the end the cell IDs range from 1 to globnIDs without gaps
        allocate(cl(nx,ny))
        CALL clustering(dat2d,globID,globID,nIDs,cl,inmissval)
        !write(*,*)"clustering: Found ",nIDs," Cells"
        !write(*,*)nIDs,globnIDs,globID
        ! periodic boundaries
        if(nIDs>0 .AND. periodic)then
          globID=globID+1-nIDs
          CALL mergeboundarycells(cl,globID,globID,nIDs,outmissval)
          !write(*,*)"perbound: Found ",nIDs," Cells"
          !write(*,*)nIDs,globnIDs,globID
        end if
        ! delete small clusters/cells
        if(nIDs>0 .AND. minarea>0)then
          globID=globID+1-nIDs
          CALL delsmallcells(cl,globID,globID,nIDs,outmissval)
          !write(*,*)"dellsmall: Found ",nIDs," Cells"
          !write(*,*)nIDs,globnIDs,globID
        end if
        if(nIDs.ne.0)globID=globID+1
        globnIDs=globnIDs+nIDs
        deallocate(dat2d)
        !write(*,*)nIDs,globnIDs,globID
    
        ! reshape array for writing to nc
        allocate(dat(nx*ny))
        CALL reshapeF2d(cl,nx,ny,dat)
        deallocate(cl)
    
        ! write time step to output file
        nmiss2=nmiss1
        CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
        deallocate(dat)
      end do
    
      ! close input and output
      CALL gridDestroy(gridID2)
      CALL vlistDestroy(vlistID2)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)
    
      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,*)"---------"
      write(*,*)"Cells found:",globnIDs
      write(*,*)"---------"
      
      if(globnIDs.eq.0)then
        write(*,*)"No cells found! We stop here!"
        write(*,*)"======================================="
        write(*,*)"======= FINISHED CELL DETECTION ======="
        write(*,*)"======================================="
        stop
      end if
    
      write(*,*)"======================================="
      write(*,*)"======= FINISHED CELL DETECTION ======="
      write(*,*)"======================================="
    
    end subroutine docelldetection
    
    subroutine clustering(data2d,startID,finID,numIDs,tcl,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres
      use ncdfpars, only : outmissval
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(out) :: finID,numIDs
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2)
      real(kind=stdfloattype), intent(in) :: data2d(nx,ny),missval
      real(kind=stdfloattype),intent(out) :: tcl(nx,ny)
      logical :: mask(nx,ny)
    
      ! initialize variables and arrays
      tcl=outmissval
      mask=.false.
      clID=startID
      numIDs=0
    
      ! mask values higher than threshold and if not missing value
      do y=1,ny
        do x=1,nx
          if(data2d(x,y)>thres .AND. data2d(x,y).ne.missval)then
          mask(x,y)=.true.
          end if
        end do
      end do
      
      ! check if there are any gridpoints to cluster
      if(ANY(mask))then
      
        ! assign IDs to continous cells
        do y=1,ny
          do x=1,nx
            neighb=outmissval
            if(mask(x,y))then
              ! gather neighbouring IDs; left,up
              if(x.ne.1)  neighb(1)=tcl((x-1),y)
              if(y.ne.1)  neighb(2)=tcl(x,(y-1))
              ! check if there is NO cluster around the current pixel; create new one
              if(ALL(neighb==outmissval))then
                tcl(x,y)=clID
                clID=clID+1
                numIDs=numIDs+1
              else
                ! both neighbours are in the same cluster
                if(neighb(1)==neighb(2).and.neighb(1).ne.outmissval)then
                  tcl(x,y)=neighb(1)
                end if
                ! both neighbors are in different clusters but none of them is (-999)
                if(neighb(1).ne.neighb(2) .and. neighb(1).ne.outmissval .and. neighb(2).ne.outmissval)then
                  numIDs=numIDs-1
                  tcl(x,y)=MINVAL(neighb)
                  ! update the existing higher cluster with the lowest neighbour
                  do cony=1,ny
                    do conx=1,nx
                      if(tcl(conx,cony)==MAXVAL(neighb))tcl(conx,cony)=MINVAL(neighb)
                    end do
                  end do
                end if
                ! both neighbors are in different clusters but ONE of them is empty(-999)
                if(neighb(1).ne.neighb(2) .and. (neighb(1)==outmissval .or. neighb(2)==outmissval))then
                  tcl(x,y)=MAXVAL(neighb)
                end if
              end if
            end if
          end do
        end do
        
        ! gather IDs and rename to gapless ascending IDs
        if(numIDs>0)then
          allocate(allIDs(numIDs))
          allIDs=outmissval
          clID=startID-1
          tp=1
          do y=1,ny
            do x=1,nx
              if(.NOT.ANY(allIDs==tcl(x,y)) .AND. tcl(x,y).ne.outmissval)then
                allIDs(tp)=tcl(x,y)
                tp=tp+1
              end if
            end do
          end do
      
          do i=1,tp-1
            clID=clID+1
            do y=1,ny
              do x=1,nx
                if(tcl(x,y)==allIDs(i))then
                  tcl(x,y)=clID
                end if
              end do
            end do
          end do
          deallocate(allIDs)
        end if
        
      end if
      ! return final cluster ID
      finID=clID
    end subroutine clustering

    subroutine mergeboundarycells(data2d,startID,finID,numIDs,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(inout) :: numIDs
      integer, intent(out) :: finID
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2)
      real(kind=stdfloattype), intent(in) :: missval
      real(kind=stdfloattype),intent(inout) :: data2d(nx,ny)
    
      ! initialize variables and arrays
    
      ! check in y direction if there are cells connected beyond boundaries
      do y=1,ny
        neighb=missval
        ! gather neighbouring IDs
        if(data2d(nx,y).ne.missval .AND. data2d(1,y).ne.missval)then
          neighb(1)=data2d(nx,y) ! this is the ID of the current gridpoint
          neighb(2)=data2d(1,y)  ! this is the ID of the neighbouring point over the edge of x
          if(neighb(1).ne.neighb(2))then
            ! now iterate the complete matrix and rename higher cluster ID to the lower ID
            do cony=1,ny
              do conx=1,nx
                if(data2d(conx,cony)==MAXVAL(neighb))data2d(conx,cony)=MINVAL(neighb)
              end do
            end do
            if(verbose)write(*,*)"y Cluster # ",MAXVAL(neighb)," replaced  by",MINVAL(neighb)
            ! one cluster was deleted
            numIDs=numIDs-1            
          end if
        end if
      end do
      ! check in x direction if there are cells connected beyond boundaries
      do x=1,nx
        neighb=missval
        ! gather neighbouring IDs
        if(data2d(x,ny).ne.missval .AND. data2d(x,1).ne.missval)then
          neighb(1)=data2d(x,ny) ! this is the ID of the current gridpoint
          neighb(2)=data2d(x,1)  ! this is the ID of the neighbouring point over the edge of y
          if(neighb(1).ne.neighb(2))then
            ! now iterate the complete matrix and rename higher cluster ID to the lower ID
            do cony=1,ny
              do conx=1,nx
                if(data2d(conx,cony)==MAXVAL(neighb))data2d(conx,cony)=MINVAL(neighb)
              end do
            end do
            if(verbose)write(*,*)"x Cluster # ",MAXVAL(neighb)," replaced  by",MINVAL(neighb)
            ! one cluster was deleted
            numIDs=numIDs-1            
          end if
        end if
      end do
      
      ! gather IDs and rename to gapless ascending IDs
      allocate(allIDs(numIDs))
      allIDs=-999
      tp=1
      do y=1,ny
        do x=1,nx
          if(.NOT.ANY(allIDs==data2d(x,y)) .AND. data2d(x,y).ne.missval)then
            allIDs(tp)=data2d(x,y)
            tp=tp+1
          end if
        end do
      end do
  
      clID=startID-1
      do i=1,tp-1
        clID=clID+1
        do y=1,ny
          do x=1,nx
            if(data2d(x,y)==allIDs(i))then
              data2d(x,y)=clID
            end if
          end do
        end do
      end do
      deallocate(allIDs)
        
      ! return final cluster ID
      finID=clID
    end subroutine mergeboundarycells
    
    subroutine delsmallcells(data2d,startID,finID,numIDs,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres,minarea
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(inout) :: numIDs
      integer, intent(out) :: finID
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2),area(numIDs)
      real(kind=stdfloattype), intent(in) :: missval
      real(kind=stdfloattype),intent(inout) :: data2d(nx,ny)
    
      ! calculate each cells area
      area=0
      do x=1,nx
        do y=1,ny
          if(data2d(x,y).ne.missval)then
          !write(*,*)data2d(x,y),x,y
          area(INT(data2d(x,y))+1-startID) = area(INT(data2d(x,y))+1-startID) + 1
          end if
        end do
      end do
      
      ! how many clusters will be deleted?
      do i=1,numIDs
        if(area(i)<minarea)numIDs=numIDs-1 
      end do
      
      ! delete clusters with area smaller than minarea
      do y=1,ny
        do x=1,nx
          if(data2d(x,y).ne.missval)then
            if(area(INT(data2d(x,y))+1-startID)<minarea)then
              data2d(x,y)=missval
              if(verbose)write(*,*)"Cluster # ",i," deleted"
            end if
          end if
        end do
      end do      
      
      if(numIDs>0)then ! otherwise all clusters were deleted!
        ! gather IDs and rename to gapless ascending IDs
        allocate(allIDs(numIDs))
        allIDs=missval
        tp=1
        do y=1,ny
          do x=1,nx
            if(.NOT.ANY(allIDs==data2d(x,y)) .AND. data2d(x,y).ne.missval)then
              allIDs(tp)=data2d(x,y)
              tp=tp+1
            end if
          end do
        end do
    
        clID=startID-1
        do i=1,tp-1
          clID=clID+1
          do y=1,ny
            do x=1,nx
              if(data2d(x,y)==allIDs(i))then
                data2d(x,y)=clID
              end if
            end do
          end do
        end do
        deallocate(allIDs)
        
        ! return final cluster ID
        finID=clID
      else
        finID=startID
      end if
    
    end subroutine delsmallcells

end module celldetection
