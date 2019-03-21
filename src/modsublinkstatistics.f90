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
module sublinkstatistics

  use globvar

  implicit none

  contains
    subroutine calcsublinkstatistics

      use globvar

      implicit none

      integer :: cltp(globnIDs)             ! counter of subcells for each cell

      write(*,*)"======================================="
      write(*,*)"========= SUBLINK STATISTICS =========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== calc number of sublinks per cell ..."
      write(*,*)"---------"

      ! now we can find out how many subcells a cell has
      allocate(clIDnsub(globnIDs))
      clIDnsub=0
      do i=1,globsubnIDs
        clIDnsub(sublinks(i))=clIDnsub(sublinks(i))+1
      end do

      write(*,*)"======================================="
      write(*,*)"=== find subcells for each cell ..."
      write(*,*)"---------"

      ! which subcells belong to a cell
      allocate(clIDsub(globnIDs,MAXVAL(clIDnsub)))
      clIDsub=-1
      cltp=1 ! this is an array with 1 dimension, see above!
      do i=1,globnIDs
        do j=1,globsubnIDs
          if(sublinks(j)==i)then
            clIDsub(i,cltp(i))=j
            cltp(i)=cltp(i)+1
          end if
        end do
      end do
      
      write(*,*)"======================================="
      write(*,*)"=== writing links to sublinks.txt ..."
      write(*,*)"---------"
      
      ! write the links to file
      open(unit=1,file="sublinks.txt",action="write",status="replace")
      write(1,*)"   subclID       clID"
      do i=1,globsubnIDs
        write(1,'(2i11)')subclIDs(i),sublinks(i)
      end do
      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"=== writing detailed subcell stats for each cell to cells_subcell_stats.txt ..."
      write(*,*)"---------"
      
      ! write detailed subcell statistics for each cell, similar to meta track statistics
      open(unit=1,file="cell_subcell_stats.txt",action="write",status="replace")
      do i=1,globnIDs
        ! write header for this cell
        write(1,'(1a4,1i12,1i4)')"### ",i,clIDnsub(i)
        write(1,*)"   subclID  subtsclID  subclarea     subclcmassX"//&
        &"     subclcmassY    subwclcmassX    subwclcmassY          subpeakVal"//&
        &"            subavVal subtouchb       date(YYYYMMDD)         time(hhmmss)"
        do k=1,clIDnsub(i)
          write(1,'(3i11,4f16.6,2f20.12,1L10,2i21.6)')subclIDs(clIDsub(i,k)),subtsclID(clIDsub(i,k)), &
          & subclarea(clIDsub(i,k)),subclcmass(clIDsub(i,k),:), &
          & subwclcmass(clIDsub(i,k),:),subclpint(clIDsub(i,k)),subclavint(clIDsub(i,k)), &
          & subtouchb(clIDsub(i,k)),subcldate(clIDsub(i,k)),subcltime(clIDsub(i,k))
        end do
      end do

      write(*,*)"======================================="
      write(*,*)"===== FINISHED SUBLINK STATISTICS ====="
      write(*,*)"======================================="

    end subroutine calcsublinkstatistics

end module sublinkstatistics
