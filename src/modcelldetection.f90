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
      real(kind=8), allocatable :: cl(:,:)
      
      ! data arrays
      real(kind=8), allocatable :: dat(:)          ! array for reading float from nc
      real(kind=8), allocatable :: dat2d(:,:)      ! array for doing the clustering
      
      globnIDs=0
      globID=1
      
      write(*,*)"======================================="
      write(*,*)"======== START CELL DETECTION ========="
      write(*,*)"======================================="
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to input file..."
    
      ! Get initial Information about grid and timesteps of both files
      CALL datainfo(ifile)
    
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
    
      ! Get information about longitudes and latitudes and levels
      allocate(xvals(0:nx-1))
      allocate(yvals(0:ny-1))
      if(verbose)write(*,*)"Arrays for grid values successfully allocated!"
      nlev=zaxisInqSize(zaxisID1)
      allocate(levels(1:nlev))
      call zaxisInqLevels(zaxisID1,levels)
      nblon=gridInqXvals(gridID1,xvals)
      nblat=gridInqYvals(gridID1,yvals)
      diflon=(xvals(nx-1)-xvals(0))/(nblon-1)
      diflat=(yvals(ny-1)-yvals(0))/(nblat-1)
      level=zaxisInqLevel(zaxisID1,levelID)
      CALL vlistInqVarUnits(vlistID1,varID1,vunit)
      CALL gridInqXunits(gridID1,xunit)
      CALL gridInqYunits(gridID1,yunit)
      inmissval=vlistInqVarMissval(vlistID1,varID1)
      call vlistInqVarName(vlistID1,varID1,vname)
      ! extract dates and times for all time steps
      allocate(vdate(ntp))
      allocate(vtime(ntp))
      do tsID=0,(ntp-1)
        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)
        ! read date and time
        vdate(tsID+1) = taxisInqVdate(taxisID1)
        vtime(tsID+1) = taxisInqVtime(taxisID1)
      end do
      
      write(*,*)"======================================="
      write(*,*)"=== INPUT SUMMARY:"
      write(*,*)"---------"
      write(*,*)"Input   : ",trim(ifile)
      write(*,*)"---------"
      write(*,'(A,1a12)')" VAR        : ",trim(vname)
      write(*,'(A,1a12)')" Unit       : ",trim(vunit)
      write(*,'(A,1f12.2)')" MissVal    : ",inmissval
      write(*,'(A,1i12)')" NX         : ",nx
      write(*,'(A,1f12.2)')" MIN X      : ",xvals(0)
      write(*,'(A,1f12.2)')" MAX X      : ",xvals(nblon-1)
      write(*,'(A,1f12.2)')" DIF X      : ",diflon
      write(*,'(A,1a12)')" Unit       : ",trim(xunit)
      write(*,'(A,1i12)')" NY         : ",ny
      write(*,'(A,1f12.2)')" MIN Y      : ",yvals(0)
      write(*,'(A,1f12.2)')" MAX Y      : ",yvals(nblat-1)
      write(*,'(A,1f12.2)')" DIF Y      : ",diflat
      write(*,'(A,1a12)')" Unit       : ",trim(yunit)
      write(*,'(A,1i8,1i8.6)')" START DATE : ",vdate(1),vtime(1)
      write(*,'(A,1i8,1I8.6)')" END DATE   : ",vdate(ntp),vtime(ntp)
      write(*,'(A,1i12)')" NTSTEPS    : ",ntp
      write(*,'(A,1i12)')" TSTEP      : ",tstep
      write(*,'(A,1i12)')" NLEV       : ",nlev
      write(*,'(A,1f12.2)')" SELLEV     : ",level
      write(*,*)"---------"
    
      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(outfile)
      write(*,*)"---------"
    
      !! open new nc file for results
      ! define grid
      gridID2=gridCreate(GRID_GENERIC, nx*ny)
      CALL gridDefXsize(gridID2,nx)
      CALL gridDefYsize(gridID2,ny)
      CALL gridDefXvals(gridID2,xvals)
      CALL gridDefYvals(gridID2,yvals)
      CALL gridDefXunits(gridID2,TRIM(xunit))
      CALL gridDefYunits(gridID2,TRIM(yunit))
      zaxisID2=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID2, level)
      ! define variables
      outmissval=-999.D0
      vlistID2=vlistCreate()
      varID2=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID2,varID2,"cellID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each cell")
      CALL vlistDefVarUnits(vlistID2,varID2,"-")
      CALL vlistDefVarMissval(vlistID2,varID2,outmissval)
      CALL vlistDefVarDatatype(vlistID2,varID2,DATATYPE_INT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite(trim(outfile),FILETYPE_NC)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)
    
      ! Get data
      write(*,*)"======================================="
      write(*,*)"=== Find continous cells ..."
      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1)then
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
        if(nIDs.ne.0)globID=globID+1
        globnIDs=globnIDs+nIDs
        deallocate(dat2d)
    
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
    
      write(*,*)"======================================="
      write(*,*)"======= FINISHED CELL DETECTION ======="
      write(*,*)"======================================="
    
    end subroutine docelldetection
    
    subroutine clustering(data2d,startID,finID,numIDs,tcl,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres,periodic
      use ncdfpars, only : outmissval
      
      implicit none
      integer, intent(in) :: startID
      integer, intent(out) :: finID,numIDs
      integer, allocatable :: allIDs(:)
      integer :: conx,cony,neighb(2)
      real(kind=8), intent(in) :: data2d(nx,ny),missval
      real(kind=8),intent(out) :: tcl(nx,ny)
      logical :: mask(nx,ny)
    
      ! initialize variables and arrays
      tcl=outmissval
      mask=.false.
      clID=startID
      numIDs=0
    
      ! mask values higher than threshold and if not missing value
      do y=1,ny
        do x=1,nx
          if(data2d(x,y)>thres .AND. data2d(x,y).ne.missval)mask(x,y)=.true.
        end do
      end do
      
      ! check if there are any gridpoints to cluster
      if(ANY(mask))then
      
        ! assign IDs to continous cells
        do y=1,ny
          do x=1,nx
            neighb=-999
            if(mask(x,y))then
              ! gather neighbouring IDs; 1= left, 2= up
              if(periodic)then
                if(x.ne.1)then
                  neighb(1)=tcl((x-1),y)
                else
                  neighb(1)=tcl(nx,y)
                end if
                if(y.ne.1)then
                  neighb(2)=tcl(x,(y-1))
                else
                  neighb(2)=tcl(x,ny)
                end if
              else
                if(x.ne.1)  neighb(1)=tcl((x-1),y)
                if(y.ne.1)  neighb(2)=tcl(x,(y-1))
              end if
              ! check if there is NO cluster around the current pixel; create new one
              if(ALL(neighb==-999))then
                tcl(x,y)=clID
                clID=clID+1
                numIDs=numIDs+1
              else
                ! both neighbours are in the same cluster
                if(neighb(1)==neighb(2).and.neighb(1).ne.-999)then
                  tcl(x,y)=neighb(1)
                end if
                ! both neighbors are in different clusters but none of them is (-999)
                if(neighb(1).ne.neighb(2) .and. neighb(1).ne.-999 .and. neighb(2).ne.-999)then
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
                if(neighb(1).ne.neighb(2) .and. (neighb(1)==-999 .or. neighb(2)==-999))then
                  tcl(x,y)=MAXVAL(neighb)
                end if
              end if
            end if
          end do
        end do
      
        ! gather IDs and rename to gapless ascending IDs
        if(numIDs>0)then
          allocate(allIDs(numIDs))
          allIDs=-999
          clID=startID-1
          tp=1
          do y=1,ny
            do x=1,nx
              if(.NOT.ANY(allIDs==tcl(x,y)) .AND. tcl(x,y).ne.-999)then
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

end module celldetection
