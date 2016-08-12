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
      
      write(*,*)"======================================="
      write(*,*)"=========== LINK STATISTICS ==========="
      write(*,*)"======================================="
    
      write(*,*)"======================================="
      write(*,*)"=== calc number of links per cell ..."
      write(*,*)"---------"
    
      ! allocate
      allocate(nbw(globnIDs),nfw(globnIDs))
      nbw=0
      nfw=0
    
      ! calc nbw/nfw
      do i=1,globnIDs
        if(tsclID(i).ne.1 .AND. iclIDloc(i).ne.-1)then
          do k=1,iclIDloc(i)
            if(links(i,k))nbw(i)=nbw(i)+1
          end do
        end if
        if(tsclID(i).ne.ntp .AND. iclIDloc(i).ne.-1)then
          do k=iclIDloc(i),maxnIDs
            if(links(i,k))nfw(i)=nfw(i)+1
          end do
        end if
      end do

      if(.NOT.advcor .OR. adviter>nadviter)then    
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
      end if
      
      write(*,*)"======================================="
      write(*,*)"====== FINISHED LINK STATISTICS ======="
      write(*,*)"======================================="
    
    end subroutine calclinkstatistics
    
end module linkstatistics
