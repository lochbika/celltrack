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
module globvar
  implicit none

  ! variables for data types
  integer, parameter :: stdclen = 1000
  integer, parameter :: stdfloattype = 8
  
  ! command line options
  character(len=stdclen) :: ifile            ! input file name
  real(kind=stdfloattype) :: thres                      ! threshold for cell detection
  real(kind=stdfloattype) :: pherevap                   ! pheromone evaporation value
  logical :: verbose                         ! verbosity
  logical :: lout                            ! switch for file output of the logical link matrix
  character(len=stdclen) :: ivar             ! name of the input variable (must match the name in NCDF file)
  integer :: levelID                         ! netcdf level ID of input variable
  integer :: nants                           ! number of ants ie agents during mainstream detection
  integer :: maxnants                        ! maximum number of ants
  integer :: nruns                           ! number of iterations during mainstream detection
  integer :: alpha,beta                      ! weighting parameters for action choice rule in ACO path algorhitm
  integer :: rseed                           ! the seed for the random number generator
  logical :: advcor                          ! switch for advection correction
  logical :: metanc                          ! switch for tracks netcdf output
  logical :: tracknc                         ! switch for meta tracks netcdf output
  logical :: nometa                          ! completely switch off meta track and mainstream routines
  logical :: periodic                        ! switch for periodic boundaries in input fields
  integer :: coarsex,coarsey                 ! factor for coarse graining of the grid for advection correction
  integer :: tstep                           ! timestep of input data in seconds
  integer :: minarea                         ! minimum area for clusters in grid points
  real(kind=stdfloattype) :: sigma                      ! the std dev for the gaussian blur before subcelldetection
  integer :: truncate                        ! truncation(span) for the gaussian blur before subcelldetection
  
  ! variables containing information about the domain
  real(kind=stdfloattype) :: level,diflon,diflat
  real(kind=stdfloattype), allocatable :: xvals(:),yvals(:)!,levels(:)
  ! Variables for general information about dimensions of input/output fields
  integer :: nx,ny,ntp,nlev,nblon,nblat
  integer, allocatable :: vdate(:),vtime(:)

  ! variables containing information about cells
  integer :: globnIDs
  integer, allocatable :: clIDs(:),tsclID(:),cldate(:),cltime(:)
  logical, allocatable :: touchb(:)
  ! Variables used for cell statistics
  integer, allocatable :: clarea(:)
  real(kind=stdfloattype), allocatable :: clpint(:),clavint(:),clcmass(:,:),wclcmass(:,:)

  ! variables containing information about *SUB*cells
  integer :: globsubnIDs
  integer, allocatable :: subclIDs(:),subtsclID(:),subcldate(:),subcltime(:)
  logical, allocatable :: subtouchb(:)
  real(kind=stdfloattype), allocatable :: kernel(:,:)
  ! Variables used for cell statistics
  integer, allocatable :: subclarea(:)
  real(kind=stdfloattype), allocatable :: subclpint(:),subclavint(:),subclcmass(:,:),subwclcmass(:,:)
  
  ! Variables for cell linking
  integer, allocatable :: nbw(:),nfw(:),clink(:),nlinks(:),links(:,:)
  integer, allocatable :: ltype(:,:) ! 0 for forward; 1 for backward link
  logical, allocatable :: tsALLna(:)
  
  ! variables for linking SUBcells to cells
  integer, allocatable :: sublinks(:) ! each subcell can only be linked to one cell :)
  integer, allocatable :: clIDnsub(:) ! how many subcells does a cell have?
  integer, allocatable :: clIDsub(:,:) ! which subcells are linked to which cell?

  ! variables for advection correction
  character(len=stdclen) :: vfile                         ! basename for velocity fields
  integer :: adviter,nadviter                         ! the current and number of iteration of advection correction
  real(kind=stdfloattype), allocatable :: uvfield2d(:,:),uvfield(:),vvfield2d(:,:),vvfield(:)   ! the velocity fields in 1d and 2d
  integer :: vnx,vny
  real(kind=stdfloattype), allocatable :: vxvals(:),vyvals(:),vclx(:),vcly(:)
  integer, allocatable :: vclxindex(:),vclyindex(:)   ! the indices of the neares gridpoints for each cell (x and y coord)
  real(kind=stdfloattype) :: maxvel ! define a maximum velocity for cells

  ! Variables used during tracking
  integer :: ntracks,ncleantr,maxtrlen
  integer, allocatable :: alltracks(:,:),trtypes(:)
  logical, allocatable :: nobounds(:)

  ! Variables used for track statistics and summary
  real(kind=stdfloattype), allocatable :: trpint(:),travint(:)
  integer, allocatable :: trdur(:),trpinttime(:)

  ! variables for meta tracking
  integer :: nmeta,maxmetalen
  integer, allocatable :: allmeta(:,:),nkinds(:,:),allcon(:,:,:)

  ! variables for meta track statistics and summary
  logical, allocatable :: mnobounds(:)
  integer, allocatable :: clmeta(:),metadur(:)
  real(kind=stdfloattype), allocatable :: metavalsum(:)
  
  ! variables for meta track mainstream detection
  integer, allocatable :: allmainstream(:,:)
  logical, allocatable :: mstrnobounds(:)

  ! variables for meta track mainstream statistics
  integer, allocatable :: clmetamstr(:),mstrdur(:)
  real(kind=stdfloattype), allocatable :: mstrvalfrac(:),mstrvalsum(:)

  ! auxiliary
  integer :: outstep,status,riostat
  character(len=stdclen) :: filename
  character(len=stdclen) :: outfile
  character(len=stdclen) :: suboutfile
  character(len=stdclen) :: blurfile
  real(kind=stdfloattype) :: pi=3.141592653589793238462643383279502884197169399373510

  ! random number
  real(kind=stdfloattype) :: rnum
  integer,allocatable :: rseeda(:)
  integer :: rsize

  ! iterators
  integer :: clID,tp,ltp,x,y,i,k,j,l,n,p,ctrack,cmeta

end module globvar
