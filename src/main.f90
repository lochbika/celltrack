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
  use celldetection
  use cellstatistics
  use celllinking
  use linkstatistics
  use buildtracks
  use trackstatistics
  use buildmetatracks
  use metatrackstatistics
  use mainstreamdetection
  use advectioncorrection

  implicit none

  ! default internal options
  outstep=50
  ! track buffer length
  maxtrlen=5000

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
  !===== START ADVECTION CORRECTION ======
  !=======================================

  if(advcor)CALL doadvectioncorrection()

  !=======================================
  !==== FINISHED ADVECTION CORRECTION ====
  !=======================================
  
  !=======================================
  !======== START CELL DETECTION =========
  !=======================================

  CALL docelldetection()

  !=======================================
  !======= FINISHED CELL DETECTION =======
  !=======================================

  !=======================================
  !=========== CELL STATISTICS ===========
  !=======================================

  CALL calccellstatistics()

  !=======================================
  !========= FINISHED STATISTICS =========
  !=======================================

  !=======================================
  !=========== STARTED LINKING ===========
  !=======================================

  CALL linking()

  !=======================================
  !========== FINISHED LINKING ===========
  !=======================================

  !=======================================
  !=========== LINK STATISTICS ===========
  !=======================================

  CALL calclinkstatistics()

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
  CALL writetracks()

  !=======================================
  !====== FINISHED TRACK STATISTICS ======
  !=======================================

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
  CALL writemetatracks()

  !=======================================
  !==== FINISHED META TRACK STATISTICS ===
  !=======================================

  !=======================================
  !======== MAINSTREAM DETECTION =========
  !=======================================

  CALL domainstreamdetection()
  CALL writemetatracksmainstream()

  !=======================================
  !==== FINISHED MAINSTREAM DETECTION ====
  !=======================================

  write(*,*)"======================================="
  write(*,*)"============ FINAL SUMMARY ============"
  write(*,*)"======================================="

  write(*,*)"---------"
  write(*,*)"Cells                  : ",globnIDs
  write(*,*)"Clean tracks           : ",ncleantr
  write(*,*)"Meta tracks            : ",nmeta
  write(*,*)"---------"
  write(*,*)"Total number of tracks : ",nmeta+ncleantr

  write(*,*)"======================================="
  write(*,*)"================= EXIT ================"
  write(*,*)"======================================="

end program celltrack
