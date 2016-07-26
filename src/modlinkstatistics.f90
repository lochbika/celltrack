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
module linkstatistics
  
  use globvar

  implicit none
  
  contains
    subroutine calclinkstatistics
      
      use globvar
     
      implicit none
      
      integer, allocatable :: linksInt(:,:)
    
      write(*,*)"======================================="
      write(*,*)"=========== LINK STATISTICS ==========="
      write(*,*)"======================================="
    
      write(*,*)"======================================="
      write(*,*)"=== calc number of links per cell ..."
      write(*,*)"---------"
    
      ! allocate
      allocate(nbw(globnIDs),nfw(globnIDs))
      allocate(linksInt(globnIDs,maxnIDs))
      nbw=0
      nfw=0
    
      ! convert logical links to integer links
      do i=1,globnIDs
        do k=1,maxnIDs
          if(links(i,k))then
            linksInt(i,k)=1
          else
            linksInt(i,k)=0
          end if
        end do
      end do
    
      ! calc nbw/nfw
      do i=1,globnIDs
        if(tsclID(i).ne.1 .AND. iclIDloc(i).ne.-1)nbw(i)=SUM(linksInt(i,1:iclIDloc(i)))
        if(tsclID(i).ne.ntp .AND. iclIDloc(i).ne.-1)nfw(i)=SUM(linksInt(i,iclIDloc(i):maxnIDs))
      end do
    
      write(*,*)"======================================="
      write(*,*)"=== write link statistics to file links_stats.txt ..."
      write(*,*)"---------"
    
      ! write this to file
      open(unit=1,file="links_stats.txt",action="write",status="replace")
      write(1,*)"     clID       nbw       nfw"
      do i=1,globnIDs
        write(1,'(3i10)')clIDs(i),nbw(i),nfw(i)
      end do
      close(unit=1)
      deallocate(linksInt)
      
      write(*,*)"======================================="
      write(*,*)"====== FINISHED LINK STATISTICS ======="
      write(*,*)"======================================="
    
    end subroutine calclinkstatistics
    
end module linkstatistics
