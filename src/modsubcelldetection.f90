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
module subcelldetection
  
  use globvar

  implicit none
  
  
  contains 
    subroutine dosubcelldetection
      
      use globvar
      use ncdfpars
      
      implicit none
      
      include 'cdi.inc'      
      
      integer :: nIDs,globID,truncate
      real(kind=8) :: sigma
      real(kind=8), allocatable :: cl(:,:)
      
      ! data arrays
      real(kind=8), allocatable :: dat(:)          ! array for reading float from nc
      real(kind=8), allocatable :: dat2d(:,:)      ! array for doing the clustering
      
      globsubnIDs=0
      globID=1
      sigma=1.5
      truncate=2
      
      write(*,*)"======================================="
      write(*,*)"====== START SUBCELL DETECTION ========"
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
    
      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(suboutfile)
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
      vlistID2=vlistCreate()
      varID2=vlistDefVar(vlistID2,gridID2,zaxisID2,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID2,varID2,"cellID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each cell")
      CALL vlistDefVarUnits(vlistID2,varID2,"-")
      CALL vlistDefVarMissval(vlistID2,varID2,outmissval)
      CALL vlistDefVarDatatype(vlistID2,varID2,CDI_DATATYPE_FLT32)
      ! copy time axis from input
      taxisID2=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID2,taxisID2)
      ! Open the dataset for writing
      streamID2=streamOpenWrite(trim(suboutfile),CDI_FILETYPE_NC4)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if
      ! set netCDF4 compression
      CALL streamDefCompType(streamID2,CDI_COMPRESS_ZIP)
      CALL streamDefCompLevel(streamID2, 6)
      ! Assign variables to dataset
      call streamDefVList(streamID2,vlistID2)
      
      ! at this stage we get the gaussian kernel for the smoothing of the input field
      CALL gaussian_kernel(sigma,truncate)
      !write(*,*)size(kernel,1),size(kernel,2)
      !do x=-3,3
      !write(*,'(7f8.6)')kernel(x,:)
      !enddo

      ! Get data
      write(*,*)"======================================="
      write(*,*)"=== Find subcells ..."
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
        if(tsALLna(tsID+1))then
          nmiss2=nmiss1
          dat=outmissval
          cycle
        end if
        
        ! reshape array
        allocate(dat2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        deallocate(dat)
    
        ! cluster the frame
        ! it is very important that at the end the cell IDs range from 1 to globnIDs without gaps
        allocate(cl(nx,ny))
        CALL blur2d(dat2d,cl,inmissval)
        !cl=dat2d
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
      write(*,*)"Subcells found:",globsubnIDs
      write(*,*)"---------"
    
      write(*,*)"======================================="
      write(*,*)"======= FINISHED CELL DETECTION ======="
      write(*,*)"======================================="
    
    end subroutine dosubcelldetection
    
    subroutine blur2d(data2d,tcl,missval)
      
      use globvar, only : y,x,i,tp,nx,ny,thres,kernel,periodic
      use ncdfpars, only : outmissval
      
      implicit none
      integer :: conx,cony,neighb(2)
      integer :: skernx,skerny
      real(kind=8), intent(in) :: data2d(nx,ny),missval
      real(kind=8),intent(out) :: tcl(nx,ny)
      real(kind=8), allocatable :: tcltmp(:,:)      ! array for extension of the domain (handling boundaries)
      logical :: mask(nx,ny)
      
      ! get the size of the kernel
      skernx=size(kernel,1)
      skerny=size(kernel,2)
    
      ! initialize variables and arrays
      tcl=outmissval
      mask=.false.
    
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
        ! now continue differently if we have periodic boundaries or not
        allocate(tcltmp( (1-skernx):(nx+skernx) , (1-skerny):(ny+skerny) ))
        tcltmp(1:nx,1:ny) = data2d

        if(periodic)then
          !copy from the other edges
          ! left side
          tcltmp((1-skernx):0,1:ny) = data2d( ((nx-(skernx-1))):nx,1:ny)
          ! right side
          tcltmp((nx+1):(nx+skernx),1:ny) = data2d( 1:skernx,1:ny)
          ! upper side
          tcltmp(1:nx,(1-skerny):0) = data2d( 1:nx,((ny-(skerny-1))):ny)
          ! lower side
          tcltmp(1:nx,(ny+1):(ny+skerny)) = data2d( 1:nx,1:skerny)
          ! upper left
          tcltmp((1-skernx):0,(1-skerny):0) = data2d( ((nx-(skernx-1))):nx,((ny-(skerny-1))):ny)
          ! upper right
          tcltmp((nx+1):(nx+skernx),(1-skerny):0) = data2d( 1:skernx,((ny-(skerny-1))):ny)
          ! lower left
          tcltmp((1-skernx):0,(ny+1):(ny+skerny)) = data2d( ((nx-(skernx-1))):nx,1:skerny)
          ! lower right
          tcltmp((nx+1):(nx+skernx),(ny+1):(ny+skerny)) = data2d( 1:skernx,1:skerny)
        end if
                
        ! blur it with respecting the boundaries
        do y=1,ny
          do x=1,nx
            tcl(x,y) = sum( tcltmp( x-((skernx-1)/2):x+((skernx-1)/2) , y-((skerny-1)/2):y+((skerny-1)/2) ) * kernel )
          end do
        end do
        
      end if
    end subroutine blur2d
    
    subroutine subclustering(data2d,startID,finID,numIDs,tcl,missval)
      
      use globvar, only : clID,y,x,i,tp,nx,ny,thres
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
    end subroutine subclustering

end module subcelldetection
