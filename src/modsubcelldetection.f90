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

      integer :: nIDs,globID
      real(kind=stdfloattype), allocatable :: cl(:,:)

      ! data arrays
      real(kind=stdfloattype), allocatable :: dat(:)          ! array for reading float from nc
      real(kind=stdfloattype), allocatable :: dat2d(:,:)      ! array for doing the clustering
      real(kind=stdfloattype), allocatable :: subcl2d(:,:)    ! array holding the subcell IDs in 2D
      real(kind=stdfloattype), allocatable :: cells(:)        ! array for reading cell IDs from nc
      real(kind=stdfloattype), allocatable :: cells2d(:,:)    ! array for holding cell IDs in 2d

      globsubnIDs=0
      globID=1

      write(*,*)"======================================="
      write(*,*)"====== START SUBCELL DETECTION ========"
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
      write(*,*)"=== Opening connection to cells file..."

      ! Open the dataset 1
      streamID3=streamOpenRead(outfile)
      if(streamID3<0)then
         write(*,*)cdiStringError(streamID3)
         stop
      end if

      ! Set the variable IDs 1
      varID3=getVarIDbyName(outfile,"cellID")
      vlistID3=streamInqVlist(streamID3)
      gridID3=vlistInqVarGrid(vlistID3,varID3)
      taxisID3=vlistInqTaxis(vlistID3)
      zaxisID3=vlistInqVarZaxis(vlistID3,varID3)

      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(suboutfile)
      write(*,*)"---------"

      !! open new nc file for results
      ! define grid
      gridID2=gridDuplicate(gridID3)
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
      CALL vlistDefVarName(vlistID2,varID2,"subcellID")
      CALL vlistDefVarLongname(vlistID2,varID2,"unique ID of each subcell")
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

      write(*,*)"======================================="
      write(*,*)"=== CREATING OUTPUT ..."
      write(*,*)"Output  :     ",trim(blurfile)
      write(*,*)"---------"

      !! open new nc file for results
      ! define grid
      gridID4=gridDuplicate(gridID3)
      !CALL gridDefXsize(gridID4,nx)
      !CALL gridDefYsize(gridID4,ny)
      !CALL gridDefXvals(gridID4,xvals)
      !CALL gridDefYvals(gridID4,yvals)
      !CALL gridDefXunits(gridID4,TRIM(xunit))
      !CALL gridDefYunits(gridID4,TRIM(yunit))
      zaxisID4=zaxisCreate(ZAXIS_GENERIC, 1)
      CALL zaxisDefLevels(zaxisID4, level)
      ! define variables
      vlistID4=vlistCreate()
      varID4=vlistDefVar(vlistID4,gridID4,zaxisID4,TIME_VARIABLE)
      CALL vlistDefVarName(vlistID4,varID4,trim(vname))
      CALL vlistDefVarLongname(vlistID4,varID4,"")
      CALL vlistDefVarUnits(vlistID4,varID4,trim(vunit))
      CALL vlistDefVarMissval(vlistID4,varID4,outmissval)
      CALL vlistDefVarDatatype(vlistID4,varID4,CDI_DATATYPE_FLT32)
      ! copy time axis from input
      taxisID4=vlistInqTaxis(vlistID1)
      call vlistDefTaxis(vlistID4,taxisID4)
      ! Open the dataset for writing
      streamID4=streamOpenWrite(trim(blurfile),CDI_FILETYPE_NC4)
      if(streamID4<0)then
         write(*,*)cdiStringError(streamID4)
         stop
      end if
      ! set netCDF4 compression
      CALL streamDefCompType(streamID4,CDI_COMPRESS_ZIP)
      CALL streamDefCompLevel(streamID4, 6)
      ! Assign variables to dataset
      call streamDefVList(streamID4,vlistID4)

      ! at this stage we get the gaussian kernel for the smoothing of the input field
      CALL gaussian_kernel(sigma,truncate)

      ! Get data
      write(*,*)"======================================="
      write(*,*)"=== Find subcells ..."
      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1 .OR. verbose)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if

        ! Allocate arrays for data storage
        allocate(dat(nx*ny))
        allocate(cells(nx*ny))

        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)

        ! Set time step for cells file
        status=streamInqTimestep(streamID3,tsID)

        ! set time step for output file
        status=streamDefTimestep(streamID2,tsID)

        ! set time step for smoothed file
        status=streamDefTimestep(streamID4,tsID)

        ! Read time step from input
        call streamReadVarSlice(streamID1,varID1,levelID,dat,nmiss1)

        ! Read time step from cells input
        call streamReadVarSlice(streamID3,varID3,0,cells,nmiss3)

        ! cycle if field contains only missing values; but write it to output
        if(nmiss3==nx*ny)then
          nmiss2=nx*ny
          dat=outmissval
          CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
          deallocate(dat,cells)
          cycle
        end if

        ! reshape arrays
        allocate(dat2d(nx,ny))
        allocate(cells2d(nx,ny))
        CALL reshapeT2d(dat,nx,ny,dat2d)
        CALL reshapeT2d(cells,nx,ny,cells2d)
        deallocate(dat)
        deallocate(cells)

        ! cluster the frame
        allocate(cl(nx,ny))
        ! let us first replace the missing value from the input file with outmissval
        nmiss4=0 ! reset missing value counter
        do x=1,nx
          do y=1,ny
            if(dat2d(x,y)==inmissval)then
              dat2d(x,y)=outmissval
              nmiss4=nmiss4+1
            end if
          end do
        end do
        ! low pass filter the array and deallocate the original one
        CALL blur2d(dat2d,cl,outmissval)
        deallocate(dat2d)
        ! write blurred field to netcdf
        allocate(dat(nx*ny))
        CALL reshapeF2d(cl,nx,ny,dat)
        CALL streamWriteVar(streamID4,varID4,dat,nmiss4)
        deallocate(dat)
        ! now cluster subcells
        allocate(subcl2d(nx,ny))
        CALL subclustering(cl,globID,globID,nIDs,subcl2d,outmissval,cells2d)
        deallocate(cl)

        if(nIDs>0)then
          globID=globID+1
          globsubnIDs=globsubnIDs+nIDs
        end if

        ! reshape array for writing to nc
        allocate(dat(nx*ny))
        CALL reshapeF2d(subcl2d,nx,ny,dat)
        deallocate(subcl2d)

        ! write time step to output file
        nmiss2=nmiss1
        CALL streamWriteVar(streamID2,varID2,dat,nmiss2)
        deallocate(dat)
        ! deallocate cells
        deallocate(cells2d)
      end do

      ! close input and output
      !CALL gridDestroy(gridID2)
      CALL vlistDestroy(vlistID2)
      !CALL gridDestroy(gridID4)
      CALL vlistDestroy(vlistID4)
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)
      CALL streamClose(streamID3)
      CALL streamClose(streamID4)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,*)"---------"
      write(*,*)"Subcells found:",globsubnIDs
      write(*,*)"---------"

      write(*,*)"======================================="
      write(*,*)"===== FINISHED SUBCELL DETECTION ======"
      write(*,*)"======================================="

    end subroutine dosubcelldetection

    subroutine blur2d(data2d,tcl,missval)

      use globvar, only : y,x,i,tp,nx,ny,thres,kernel,periodic

      implicit none
      integer :: conx,cony,neighb(2)
      integer :: skernx,skerny
      integer :: kx,ky ! iterators over dimensions of filter window
      real(kind=stdfloattype) :: fltav ! average value inside filter window; NA values are replaced by this variable
      real(kind=stdfloattype), intent(in) :: data2d(nx,ny),missval
      real(kind=stdfloattype),intent(out) :: tcl(nx,ny)
      real(kind=stdfloattype), allocatable :: tcltmp(:,:)      ! array for extension of the domain (handling boundaries)
      real(kind=stdfloattype), allocatable :: flttmp(:,:)      ! array for temporary storing the values for the current filter region
      logical :: mask(nx,ny)

      ! get the size of the kernel
      skernx=size(kernel,1)
      skerny=size(kernel,2)

      ! initialize variables and arrays
      tcl=missval
      mask=.false.

      ! mask values higher than threshold and if not missing value
      do y=1,ny
        do x=1,nx
          if(data2d(x,y)>thres .AND. data2d(x,y).ne.missval)then
            mask(x,y)=.true.
          end if
        end do
      end do

      ! check if there are any gridpoints masked
      if(ANY(mask))then
        ! now continue differently if we have periodic boundaries or not
        allocate(tcltmp( (1-skernx):(nx+skernx) , (1-skerny):(ny+skerny) ))
        tcltmp=missval ! initially fill with NA
        tcltmp(1:nx,1:ny) = data2d ! fill in values we have from data2d

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
        ! this is the temporary array which holds the values of the filter window
        allocate(flttmp(skernx,skerny))
        do y=1,ny
          do x=1,nx
            if(.NOT.mask(x,y))cycle ! skip if this grid point is a missing value
            ! copy window to temporary array
            flttmp = tcltmp( x-((skernx-1)/2):x+((skernx-1)/2) , y-((skerny-1)/2):y+((skerny-1)/2) )
            ! calculate the average in the filter window
            fltav=0.D0
            tp=0
            do kx=1,skernx
              do ky=1,skerny
                if(flttmp(kx,ky).ne.missval)then
                  fltav=fltav+flttmp(kx,ky) ! sum it up
                  tp=tp+1
                end if
              end do
            end do
            fltav=fltav/tp ! average
            ! replace missing values with the average
            do kx=1,skernx
              do ky=1,skerny
                if(flttmp(kx,ky)==missval)flttmp(kx,ky)=fltav
              end do
            end do
            ! apply the filter to the current window and save it to the output array
            tcl(x,y) = sum( flttmp * kernel )
          end do
        end do

        ! set grid points that were not originally covered to missval
        do y=1,ny
          do x=1,nx
            if(.NOT.mask(x,y))tcl(x,y) = missval
          end do
        end do

      end if
    end subroutine blur2d

    subroutine subclustering(data2d,startID,finID,numIDs,tcl,missval,cIDs)

      use globvar, only : y,x,i,tp,nx,ny,thres,periodic

      implicit none
      integer, intent(in) :: startID
      integer, intent(out) :: finID,numIDs
      integer, allocatable :: ipartcoor(:,:),tpartcoor(:,:),partcoor(:,:),locmaxcoor(:,:)
      integer :: npart,nlocmax
      real(kind=stdfloattype), intent(in) :: data2d(nx,ny),cIDs(nx,ny),missval
      real(kind=stdfloattype),intent(out) :: tcl(nx,ny)
      real(kind=stdfloattype) :: lmaxt,neighb(9)
      integer :: dir(nx,ny)
      logical :: mask(nx,ny),locmax(nx,ny)

      ! initialize variables and arrays
      tcl=missval
      mask=.false.
      locmax=.false.
      dir=0
      numIDs=0
      npart=0
      nlocmax=0

      ! mask values if not missing value
      do y=1,ny
        do x=1,nx
          if(data2d(x,y).ne.missval)then
            mask(x,y)=.true.
          end if
        end do
      end do

      ! check if there are any gridpoints to cluster
      if(ANY(mask))then

        ! find local maxima and directions
        do y=1,ny
          do x=1,nx
            if(mask(x,y))then
              ! extract the 8 neighboring grid points
              ! IMPORTANT: we only extract values of neighbors which belong to the same cell!
              neighb=missval
              ! the center gridpoint is always assigned the same way
              neighb(5)=data2d(x,y)
              ! for the rest, we have to take care of the boundaries
              ! we are in the upper left corner
              if(x.eq.1 .AND. y.eq.1)then
                if(cIDs(x+1,y).eq.cIDs(x,y))neighb(6)=data2d(x+1,y)
                if(cIDs(x,y+1).eq.cIDs(x,y))neighb(8)=data2d(x,y+1)
                if(cIDs(x+1,y+1).eq.cIDs(x,y))neighb(9)=data2d(x+1,y+1)
                if(periodic)then
                  if(cIDs(nx,ny).eq.cIDs(x,y))neighb(1)=data2d(nx,ny)
                  if(cIDs(x,ny).eq.cIDs(x,y))neighb(2)=data2d(x,ny)
                  if(cIDs(x+1,ny).eq.cIDs(x,y))neighb(3)=data2d(x+1,ny)
                  if(cIDs(nx,y).eq.cIDs(x,y))neighb(4)=data2d(nx,y)
                  if(cIDs(nx,y+1).eq.cIDs(x,y))neighb(7)=data2d(nx,y+1)
                end if
              ! we are in the upper right corner
              else if(x.eq.nx .AND. y.eq.1)then
                if(cIDs(x-1,y).eq.cIDs(x,y))neighb(4)=data2d(x-1,y)
                if(cIDs(x-1,y+1).eq.cIDs(x,y))neighb(7)=data2d(x-1,y+1)
                if(cIDs(x,y+1).eq.cIDs(x,y))neighb(8)=data2d(x,y+1)
                if(periodic)then
                  if(cIDs(x-1,ny).eq.cIDs(x,y))neighb(1)=data2d(x-1,ny)
                  if(cIDs(x,ny).eq.cIDs(x,y))neighb(2)=data2d(x,ny)
                  if(cIDs(1,ny).eq.cIDs(x,y))neighb(3)=data2d(1,ny)
                  if(cIDs(1,y).eq.cIDs(x,y))neighb(6)=data2d(1,y)
                  if(cIDs(1,y+1).eq.cIDs(x,y))neighb(9)=data2d(1,y+1)
                end if
              ! we are in the lower left corner
              else if(x.eq.1 .AND. y.eq.ny)then
                if(cIDs(x,y-1).eq.cIDs(x,y))neighb(2)=data2d(x,y-1)
                if(cIDs(x+1,y-1).eq.cIDs(x,y))neighb(3)=data2d(x+1,y-1)
                if(cIDs(x+1,y).eq.cIDs(x,y))neighb(6)=data2d(x+1,y)
                if(periodic)then
                  if(cIDs(nx,y-1).eq.cIDs(x,y))neighb(1)=data2d(nx,y-1)
                  if(cIDs(nx,y).eq.cIDs(x,y))neighb(4)=data2d(nx,y)
                  if(cIDs(nx,1).eq.cIDs(x,y))neighb(7)=data2d(nx,1)
                  if(cIDs(x,1).eq.cIDs(x,y))neighb(8)=data2d(x,1)
                  if(cIDs(x+1,1).eq.cIDs(x,y))neighb(9)=data2d(x+1,1)
                end if
              ! we are in the lower right corner
              else if(x.eq.nx .AND. y.eq.ny)then
                if(cIDs(x-1,y-1).eq.cIDs(x,y))neighb(1)=data2d(x-1,y-1)
                if(cIDs(x,y-1).eq.cIDs(x,y))neighb(2)=data2d(x,y-1)
                if(cIDs(x-1,y).eq.cIDs(x,y))neighb(4)=data2d(x-1,y)
                if(periodic)then
                  if(cIDs(1,y-1).eq.cIDs(x,y))neighb(3)=data2d(1,y-1)
                  if(cIDs(1,y).eq.cIDs(x,y))neighb(6)=data2d(1,y)
                  if(cIDs(x-1,1).eq.cIDs(x,y))neighb(7)=data2d(x-1,1)
                  if(cIDs(x,1).eq.cIDs(x,y))neighb(8)=data2d(x,1)
                  if(cIDs(1,1).eq.cIDs(x,y))neighb(9)=data2d(1,1)
                end if
              ! we are in between upper left and lower left corner
              else if(x.eq.1 .AND. y.ne.1 .AND. y.ne.ny)then
                if(cIDs(x,y-1).eq.cIDs(x,y))neighb(2)=data2d(x,y-1)
                if(cIDs(x+1,y-1).eq.cIDs(x,y))neighb(3)=data2d(x+1,y-1)
                if(cIDs(x+1,y).eq.cIDs(x,y))neighb(6)=data2d(x+1,y)
                if(cIDs(x,y+1).eq.cIDs(x,y))neighb(8)=data2d(x,y+1)
                if(cIDs(x+1,y+1).eq.cIDs(x,y))neighb(9)=data2d(x+1,y+1)
                if(periodic)then
                  if(cIDs(nx,y-1).eq.cIDs(x,y))neighb(1)=data2d(nx,y-1)
                  if(cIDs(nx,y).eq.cIDs(x,y))neighb(4)=data2d(nx,y)
                  if(cIDs(nx,y+1).eq.cIDs(x,y))neighb(7)=data2d(nx,y+1)
                end if
              ! we are in between upper left and upper right corner
              else if(y.eq.1 .AND. x.ne.1 .AND. x.ne.nx)then
                if(cIDs(x-1,y).eq.cIDs(x,y))neighb(4)=data2d(x-1,y)
                if(cIDs(x+1,y).eq.cIDs(x,y))neighb(6)=data2d(x+1,y)
                if(cIDs(x-1,y+1).eq.cIDs(x,y))neighb(7)=data2d(x-1,y+1)
                if(cIDs(x,y+1).eq.cIDs(x,y))neighb(8)=data2d(x,y+1)
                if(cIDs(x+1,y+1).eq.cIDs(x,y))neighb(9)=data2d(x+1,y+1)
                if(periodic)then
                  if(cIDs(x-1,ny).eq.cIDs(x,y))neighb(1)=data2d(x-1,ny)
                  if(cIDs(x,ny).eq.cIDs(x,y))neighb(2)=data2d(x,ny)
                  if(cIDs(x+1,ny).eq.cIDs(x,y))neighb(3)=data2d(x+1,ny)
                end if
              ! we are in between lower left and lower right corner
              else if(y.eq.ny .AND. x.ne.1 .AND. x.ne.nx)then
                if(cIDs(x-1,y-1).eq.cIDs(x,y))neighb(1)=data2d(x-1,y-1)
                if(cIDs(x,y-1).eq.cIDs(x,y))neighb(2)=data2d(x,y-1)
                if(cIDs(x+1,y-1).eq.cIDs(x,y))neighb(3)=data2d(x+1,y-1)
                if(cIDs(x-1,y).eq.cIDs(x,y))neighb(4)=data2d(x-1,y)
                if(cIDs(x+1,y).eq.cIDs(x,y))neighb(6)=data2d(x+1,y)
                if(periodic)then
                  if(cIDs(x-1,1).eq.cIDs(x,y))neighb(7)=data2d(x-1,1)
                  if(cIDs(x,1).eq.cIDs(x,y))neighb(8)=data2d(x,1)
                  if(cIDs(x+1,1).eq.cIDs(x,y))neighb(9)=data2d(x+1,1)
                end if
              ! we are in between upper right and lower right corner
              else if(x.eq.nx .AND. y.ne.1 .AND. y.ne.ny)then
                if(cIDs(x-1,y-1).eq.cIDs(x,y))neighb(1)=data2d(x-1,y-1)
                if(cIDs(x,y-1).eq.cIDs(x,y))neighb(2)=data2d(x,y-1)
                if(cIDs(x-1,y).eq.cIDs(x,y))neighb(4)=data2d(x-1,y)
                if(cIDs(x-1,y+1).eq.cIDs(x,y))neighb(7)=data2d(x-1,y+1)
                if(cIDs(x,y+1).eq.cIDs(x,y))neighb(8)=data2d(x,y+1)
                if(periodic)then
                  if(cIDs(1,y-1).eq.cIDs(x,y))neighb(3)=data2d(1,y-1)
                  if(cIDs(1,y).eq.cIDs(x,y))neighb(6)=data2d(1,y)
                  if(cIDs(1,y+1).eq.cIDs(x,y))neighb(9)=data2d(1,y+1)
                end if
              ! NOW, this is the case when we don't have to care about the boundaries
              else
                if(cIDs(x-1,y-1).eq.cIDs(x,y))neighb(1)=data2d(x-1,y-1)
                if(cIDs(x,y-1).eq.cIDs(x,y))neighb(2)=data2d(x,y-1)
                if(cIDs(x+1,y-1).eq.cIDs(x,y))neighb(3)=data2d(x+1,y-1)
                if(cIDs(x-1,y).eq.cIDs(x,y))neighb(4)=data2d(x-1,y)
                if(cIDs(x+1,y).eq.cIDs(x,y))neighb(6)=data2d(x+1,y)
                if(cIDs(x-1,y+1).eq.cIDs(x,y))neighb(7)=data2d(x-1,y+1)
                if(cIDs(x,y+1).eq.cIDs(x,y))neighb(8)=data2d(x,y+1)
                if(cIDs(x+1,y+1).eq.cIDs(x,y))neighb(9)=data2d(x+1,y+1)
              end if

              ! so, now we know the neighboring grid points
              ! is this grid point a local maximum??
              lmaxt=-1.0
              tp=0
              do i=1,9
                if(neighb(i)>lmaxt .AND. neighb(i).ne.missval)then
                  tp=i
                  lmaxt=neighb(i)
                end if
              end do
              if(tp.eq.5)then
                locmax(x,y)=.true.
              end if
              ! save the direction
              dir(x,y)=tp
            end if
          end do
        end do

        ! we have the directions for each pixel
        ! and know which of them are local maxima
        ! lets seed particles and let them migrate to them

        ! how many local maxima and cell grid points
        do y=1,ny
          do x=1,nx
            if(mask(x,y))then
              npart=npart+1
            end if
            if(locmax(x,y))then
              nlocmax=nlocmax+1
            end if
          end do
        end do

        ! get their coordinates
        allocate(ipartcoor(npart,2),partcoor(npart,2),locmaxcoor(nlocmax,2))
        tp=0 ! counter for locmax
        i=0  ! counter for particles
        do y=1,ny
          do x=1,nx
            if(mask(x,y))then
              i=i+1
              ipartcoor(i,1)=x
              ipartcoor(i,2)=y
              partcoor(i,1)=x
              partcoor(i,2)=y
            end if
            if(locmax(x,y))then
              tp=tp+1
              locmaxcoor(tp,1)=x
              locmaxcoor(tp,2)=y
            end if
          end do
        end do

        ! now start the iterative process to move particles
        allocate(tpartcoor(npart,2))

        do i=1,1000
          tpartcoor=partcoor

          do tp=1,npart
            ! is this particle already at a local maximum? Then, don't move it!
            if(locmax( tpartcoor(tp,1),tpartcoor(tp,2) ))then
              cycle
            end if
            ! otherwise, get direction and move it!
            if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.1)then
              tpartcoor(tp,1)=tpartcoor(tp,1)-1
              tpartcoor(tp,2)=tpartcoor(tp,2)-1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.2)then
              tpartcoor(tp,2)=tpartcoor(tp,2)-1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.3)then
              tpartcoor(tp,1)=tpartcoor(tp,1)+1
              tpartcoor(tp,2)=tpartcoor(tp,2)-1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.4)then
              tpartcoor(tp,1)=tpartcoor(tp,1)-1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.6)then
              tpartcoor(tp,1)=tpartcoor(tp,1)+1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.7)then
              tpartcoor(tp,1)=tpartcoor(tp,1)-1
              tpartcoor(tp,2)=tpartcoor(tp,2)+1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.8)then
              tpartcoor(tp,2)=tpartcoor(tp,2)+1
            else if(dir( tpartcoor(tp,1),tpartcoor(tp,2) ).eq.9)then
              tpartcoor(tp,1)=tpartcoor(tp,1)+1
              tpartcoor(tp,2)=tpartcoor(tp,2)+1
            end if
            ! check if we moved across boundaries
            if(periodic)then
              if(tpartcoor(tp,1)<1)tpartcoor(tp,1)=nx
              if(tpartcoor(tp,1)>nx)tpartcoor(tp,1)=1
              if(tpartcoor(tp,2)<1)tpartcoor(tp,2)=ny
              if(tpartcoor(tp,2)>ny)tpartcoor(tp,2)=1
            else ! reset
              if(tpartcoor(tp,1)<1)tpartcoor(tp,1)=partcoor(tp,1)
              if(tpartcoor(tp,1)>nx)tpartcoor(tp,1)=partcoor(tp,1)
              if(tpartcoor(tp,2)<1)tpartcoor(tp,2)=partcoor(tp,2)
              if(tpartcoor(tp,2)>ny)tpartcoor(tp,2)=partcoor(tp,2)
            end if
          end do

          ! check if we moved any particle
          if(ALL(partcoor==tpartcoor))then
            if(verbose)write(*,*)"We finished after",i," iterations!"
            exit
          end if
          ! update partcoor
          partcoor=tpartcoor
        end do

        ! now find which particles ended up in which local maximum
        ! and assign IDs to gridpoints
        do i=1,npart
          do tp=1,nlocmax
            ! if coordinates match assign ID to gridpoint where the particle was originally
            if(partcoor(i,1).eq.locmaxcoor(tp,1) .AND. partcoor(i,2).eq.locmaxcoor(tp,2))then
              tcl( ipartcoor(i,1),ipartcoor(i,2) ) = (startID-1)+tp
              exit
            end if
            if(tp==nlocmax)write(*,*)"WARNING: A particle didn't move to a local maximum!!!"
          end do
        end do

      end if
      ! return final cluster ID and set numIDs to return to main program
      if(nlocmax>0)then
        finID=(startID-1)+nlocmax
      end if
      numIDs=nlocmax

      ! checks
      do y=1,ny
        do x=1,nx
          ! is there a subcell but no cell?
          if(tcl(x,y).ne.missval .AND. cIDs(x,y)==missval)then
            write(*,*)"WARNING: we have a SUBcell where there is no cell",cIDs(x,y),tcl(x,y)
          end if
          ! is there a cell but no subcell?
          if(tcl(x,y)==missval .AND. cIDs(x,y).ne.missval)then
            write(*,*)"WARNING: we have a cell where there is no SUBcell",cIDs(x,y),tcl(x,y)
          end if
        end do
      end do
      ! check if any subcell is part of more than one cell
      do i=startID,finID
        tp=0
        k=0
        do y=1,ny
          do x=1,nx
            if(tcl(x,y)==i .AND. cIDs(x,y).ne.missval)then
              if(tp.ne.cIDs(x,y))then
                tp=cIDs(x,y)
                k=k+1
              end if
              if(k>1)write(*,*)"SUBcell",i," is part of",k," cells",cIDs(x,y)
            end if
          end do
        end do
      end do
    end subroutine subclustering

end module subcelldetection
