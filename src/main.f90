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
  use advstats
  use celllinking
  use linkstatistics
  use buildtracks
  use trackstatistics
  use buildmetatracks
  use metatrackstatistics
  use mainstreamdetection
  use advectioncorrection
  use mainstreamstatistics
  use buffering
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
  ! name for the buffer file
  bffile="cells_buffer.nc"
  ! name for the buffered cells (with IDs) file
  bfclfile="cells_bufferID.nc"

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
  !======== START CELL DETECTION =========
  !=======================================

  CALL docelldetection()

  !=======================================
  !======= FINISHED CELL DETECTION =======
  !=======================================

  !=======================================
  !====== START BUFFER CLUSTERING ========
  !=======================================

  if(buffer>0)CALL dobuffercluster()

  !=======================================
  !===== FINISHED BUFFER CLUSTERING ======
  !=======================================

  !=======================================
  !=========== CELL STATISTICS ===========
  !=======================================

  CALL calccellstatistics()
  
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
