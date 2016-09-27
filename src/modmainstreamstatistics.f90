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
module mainstreamstatistics

  use globvar

  implicit none

  contains
    subroutine calcmainstreamstatistics

      use globvar

      implicit none

      write(*,*)"======================================="
      write(*,*)"======== MAINSTREAM STATISTICS ========"
      write(*,*)"======================================="

      ! allocate the arrays
      allocate(mstrvalfrac(nmeta),mstrvalsum(nmeta),metavalsum(nmeta),mstrdur(nmeta))
      mstrvalfrac = 0
      mstrvalsum  = 0
      metavalsum  = 0
      mstrdur     = 0

      write(*,*)"======================================="
      write(*,*)"=== calculating summary statistics ..."
      write(*,*)"---------"
      
      ! calculate the fraction of the key variable (e.g. precip) the mainstream has
      ! at first the sum for each meta track and the mainstream
      do i=1,globnIDs
        if(clmetamstr(i).ne.-1)mstrvalsum(clmetamstr(i)) = mstrvalsum(clmetamstr(i)) + ( clavint(i)*clarea(i) )
        if(clmeta(i).ne.-1)metavalsum(clmeta(i)) = metavalsum(clmeta(i)) + ( clavint(i)*clarea(i) )
      end do

      ! now the fraction
      do i=1,nmeta
        mstrvalfrac(i) = mstrvalsum(i)/metavalsum(i)
      end do
      
      ! the number of cells in a mainstream is the duration
      do i=1,globnIDs
        if(clmetamstr(i).ne.-1)mstrdur(clmetamstr(i)) = mstrdur(clmetamstr(i)) + 1
      end do

      write(*,*)"======================================="
      write(*,*)"=== write to meta_mainstream_summary.txt ..."
      write(*,*)"---------"

      ! write the stats to file
      open(unit=1,file="meta_mainstream_summary.txt",action="write",status="replace")
      write(1,*)"     metaID         dur      valSum     valFrac"
      do i=1,nmeta
        write(1,'(2i12,2f12.6)')i,mstrdur(i),mstrvalsum(i),mstrvalfrac(i)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"==== FINISHED MAINSTREAM STATISTICS ==="
      write(*,*)"======================================="

    end subroutine calcmainstreamstatistics

end module mainstreamstatistics
