!-------------------------------------------------------------------------------------------
!  ######  ######## ##       ##       ######## ########     ###     ######  ##    ##
! ##    ## ##       ##       ##          ##    ##     ##   ## ##   ##    ## ##   ##
! ##       ##       ##       ##          ##    ##     ##  ##   ##  ##       ##  ##
! ##       ######   ##       ##          ##    ########  ##     ## ##       #####
! ##       ##       ##       ##          ##    ##   ##   ######### ##       ##  ##
! ##    ## ##       ##       ##          ##    ##    ##  ##     ## ##    ## ##   ##
!  ######  ######## ######## ########    ##    ##     ## ##     ##  ######  ##    ##
!-------------------------------------------------------------------------------------------
! This is a cell tracking algorithm inspired by Moseley (2013) doi:10.1002/2013JD020868
! Copyright: Kai Lochbihler (kai.lochbihler@knmi.nl)
!
! Any help will be appreciated :)
!
program celltrack

  ! global variables
  use globvar
  ! modules
  use inputinfo
  use celldetection
  use subcelldetection
  use cellstatistics
  use subcellstatistics
  use advstats
  use celllinking
  use subcelllinking
  use linkstatistics
  use sublinkstatistics
  use buildtracks
  use trackstatistics
  use buildmetatracks
  use metatrackstatistics
  use mainstreamdetection
  use advectioncorrection
  use mainstreamstatistics
  use cellshape

  implicit none

  ! default internal options
  outstep=50
  ! track buffer length
  maxtrlen=5000
  ! current advection iteration 
  adviter=0
  ! name for the cells file
  outfile="cells.nc"  
  ! name for the *SUB*cells file
  suboutfile="subcells.nc"
  ! name for the smoothed input file
  blurfile="input_smoothed.nc"

  ! init the random seed if not specified with cli option
  if(rseed.eq.-1)then
    CALL init_random_seed()
  else
    CALL random_seed(size=rsize)
    allocate(rseeda(rsize))
    rseeda=rseed
    CALL random_seed(put=rseeda)
    deallocate(rseeda)
  end if

  !=======================================
  !============ CLI ARGUMENTS ============
  !=======================================
  
  CALL cliarguments()

  !=======================================
  !======= FINISHED CLI ARGUMENTS ========
  !=======================================

  !=======================================
  !========= INPUT INFORMATION ===========
  !=======================================
  
  CALL gatherInputInfo()
  
  !=======================================
  !======= END INPUT INFORMATION =========
  !=======================================
  
  !=======================================
  !======== START CELL DETECTION =========
  !=======================================

  CALL docelldetection()
  CALL dosubcelldetection()

  !=======================================
  !======= FINISHED CELL DETECTION =======
  !=======================================

  !=======================================
  !=========== CELL STATISTICS ===========
  !=======================================

  CALL calccellstatistics()
  CALL calcsubcellstatistics()
  
  CALL calccellshape()

  !=======================================
  !========= FINISHED STATISTICS =========
  !=======================================

  !=======================================
  !======== MORE  CELL STATISTICS ========
  !=======================================

  CALL calccellpercentiles()

  !=======================================
  !====== FINISHED MORE STATISTICS =======
  !=======================================

  !=======================================
  !===== START ADVECTION CORRECTION ======
  !=======================================

  if(advcor)CALL doadvectioncorrection()

  !=======================================
  !==== FINISHED ADVECTION CORRECTION ====
  !=======================================

  !=======================================
  !=========== STARTED LINKING ===========
  !=======================================

  CALL linking()
  CALL dosubcelllinking()

  !=======================================
  !========== FINISHED LINKING ===========
  !=======================================

  !=======================================
  !=========== LINK STATISTICS ===========
  !=======================================

  CALL calclinkstatistics()
  CALL calcsublinkstatistics()

  !=======================================
  !====== FINISHED LINK STATISTICS =======
  !=======================================

  !=======================================
  !=========== BUILDING TRACKS ===========
  !=======================================

  CALL dobuildtracks()

  !=======================================
  !====== FINISHED BUILDING TRACKS =======
  !=======================================

  !=======================================
  !========== TRACK STATISTICS ===========
  !=======================================

  CALL calctrackstatistics()
  if(tracknc)CALL writetracks()

  !=======================================
  !====== FINISHED TRACK STATISTICS ======
  !=======================================

  if(.NOT.nometa)then
    !=======================================
    !======== BUILDING META TRACKS =========
    !=======================================
  
    CALL dobuildmetatracks()
  
    !=======================================
    !==== FINISHED BUILDING META TRACKS ====
    !=======================================
  
    !=======================================
    !======== META TRACK STATISTICS ========
    !=======================================
  
    CALL calcmetatrackstatistics()
    if(metanc)CALL writemetatracks()
  
    !=======================================
    !==== FINISHED META TRACK STATISTICS ===
    !=======================================
  end if
  
  if(.NOT.nometa .AND. .NOT.nometamstr)then
    !=======================================
    !======== MAINSTREAM DETECTION =========
    !=======================================
  
    CALL domainstreamdetection()
    if(metanc)CALL writemetatracksmainstream()
  
    !=======================================
    !==== FINISHED MAINSTREAM DETECTION ====
    !=======================================
  
    !=======================================
    !======== MAINSTREAM STATISTICS ========
    !=======================================
  
    CALL calcmainstreamstatistics()
  
    !=======================================
    !==== FINISHED MAINSTREAM STATISTICS ===
    !=======================================
  end if

  write(*,*)"======================================="
  write(*,*)"============ FINAL SUMMARY ============"
  write(*,*)"======================================="

  write(*,*)"---------"
  write(*,*)"Cells                  : ",globnIDs
  write(*,*)"SUBcells               : ",globsubnIDs
  write(*,*)"Clean tracks           : ",ncleantr
  if(.NOT.nometa)write(*,*)"Meta tracks            : ",nmeta
  write(*,*)"---------"
  if(.NOT.nometa)then
    write(*,*)"Total number of tracks : ",nmeta+ncleantr
  else
    write(*,*)"Total number of tracks : ",ncleantr
  end if

  write(*,*)"======================================="
  write(*,*)"================= EXIT ================"
  write(*,*)"======================================="

end program celltrack
