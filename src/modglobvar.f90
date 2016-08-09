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

  ! command line options
  character(len=800) :: ifile,outfile        ! input file name
  real(kind=8) :: thres                      ! threshold for cell detection
  real(kind=8) :: pherevap                   ! pheromone evaporation value
  logical :: verbose                         ! verbosity
  logical :: lout                            ! switch for file output of the logical link matrix
  integer :: ivar                            ! netcdf ID of input variable
  integer :: levelID                         ! netcdf level ID of input variable
  integer :: nants                           ! number of ant ie agents during mainstream detection
  integer :: nruns                           ! number of iterations during mainstream detection
  integer :: alpha,beta                      ! weighting parameters for action choice rule in ACO path algorhitm
  integer :: rseed                           ! the seed for the random number generator
  logical :: advcor                          ! switch for advection correction

  ! variables containing information about the domain
  real(kind=8) :: level,diflon,diflat
  real(kind=8), allocatable :: xvals(:),yvals(:),levels(:)
  ! Variables for general information about dimensions of input/output fields
  integer :: nx,ny,ntp,nlev,nblon,nblat

  ! variables containing information about cells
  integer :: globnIDs
  integer, allocatable :: clIDs(:),tsclID(:)
  logical, allocatable :: touchb(:)
  ! Variables used for cell statistics
  integer, allocatable :: clarea(:),wclarea(:)
  real(kind=8), allocatable :: clpint(:),clavint(:),clcmass(:,:),wclcmass(:,:)

  ! Variables for cell linking
  integer, allocatable :: nbw(:),nfw(:),minclIDloc(:),iclIDloc(:)
  logical, allocatable :: links(:,:)
  integer :: maxnIDs

  ! variables for advection correction
  character(len=800) :: vfile                ! filename for velocity fields
  integer :: adviter,nadviter                ! the current and number of iteration of advection correction

  ! Variables used during tracking
  integer :: ntracks,ncleantr,maxtrlen
  integer, allocatable :: alltracks(:,:),trtypes(:)
  logical, allocatable :: nobounds(:)

  ! Variables used for track statistics and summary
  real(kind=8), allocatable :: trpint(:),travint(:)
  integer, allocatable :: trdur(:),trpinttime(:)

  ! variables for meta tracking
  integer :: nmeta,maxmetalen
  integer, allocatable :: allmeta(:,:),nkinds(:,:),allcon(:,:,:)

  ! variables for track statistics and summary
  logical, allocatable :: mnobounds(:)
  integer, allocatable :: clmeta(:),metadur(:),clmetamstr(:)

  ! variables for meta track mainstream detection
  integer, allocatable :: allmainstream(:,:)

  ! auxiliary
  integer :: outstep,status,riostat

  ! random number
  real(kind=8) :: rnum
  integer,allocatable :: rseeda(:)
  integer :: rsize

  ! iterators
  integer :: clID,tp,ltp,x,y,i,k,j,l,n,p,ctrack,cmeta

end module globvar
