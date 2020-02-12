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
subroutine cliarguments

  use globvar

  implicit none

  character(len=stdclen) :: argc,command
  integer :: narg,arg

  call getarg(0,command)
  narg=IARGC()
  if(narg==0)call help(trim(command))
  arg=0

  ! defaults
  verbose=.false.
  thres=0.D0
  ifile=""
  ivar=""
  levelID=0
  nants=-1
  maxnants=-1
  nruns=300
  lout=.false.
  pherevap=0.5D0
  alpha=1
  beta=2
  rseed=-1
  advcor=.false.
  nadviter=6
  coarsex=10
  coarsey=10
  tstep=-1
  metanc=.false.
  tracknc=.false.
  nometa=.false.
  nometamstr=.false.
  periodic=.false.
  maxvel=50
  minarea=0
  sigma=2
  truncate=4
  subc=.false.
  thresdir=1

  do while (arg < narg)
    arg=arg+1
    call getarg(arg,argc)

    select case (trim(argc))
    case ("-h")
      call help("")
    case ("-v")
      verbose=.true.
    case ("-lout")
      lout=.true.
    case ("-i")
      arg=arg+1
      call getarg(arg,ifile)
    case ("-thres")
      arg=arg+1
      call getarg(arg,argc)
      do i=1,stdclen
        if(argc(i:i)=="+")then
          thresdir=1
          k=i
          exit
        else if(argc(i:i)=="-")then
          thresdir=0
          k=i
          exit
        else
          thresdir=1
        endif  
      end do
      argc(k:k)=" "
      read(argc,*)thres
    case ("-var")
      arg=arg+1
      call getarg(arg,ivar)
    case ("-lev")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)levelID
    case ("-nants")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)nants
    case ("-maxnants")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)maxnants
    case ("-nruns")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)nruns
      if(nruns<5)nruns=5
    case ("-rho")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)pherevap
    case ("-rseed")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)rseed
    case ("-alpha")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)alpha
    case ("-beta")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)beta
    case ("-nadviter")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)nadviter
    case ("-cx")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)coarsex
    case ("-cy")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)coarsey
    case ("-tstep")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)tstep
    case ("-advcor")
      advcor=.true.
    case ("-subc")
      subc=.true.
    case ("-tracknc")
      tracknc=.true.
    case ("-metanc")
      metanc=.true.
    case ("-nometa")
      nometa=.true.
    case ("-nometamstr")
      nometamstr=.true.
    case ("-perbound")
      periodic=.true.
    case ("-maxv")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)maxvel
    case ("-minarea")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)minarea
    case ("-sigma")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)sigma
    case ("-trunc")
      arg=arg+1
      call getarg(arg,argc)
      read(argc,*)truncate
    case DEFAULT
      call help(trim(command)//": ERROR: unknown argument: "//trim(argc))
    end select

    cycle
  enddo

  ! check if all necessary variables are set
  if(trim(ifile)=="")call help("No input file selected!")
  if(tstep==-1 .AND. advcor)call help("Missing argument tstep... necessary if you want to do advection correction")

end subroutine cliarguments
