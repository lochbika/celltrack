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
module buildtracks

  use globvar

  implicit none

  contains
    subroutine dobuildtracks

      use globvar

      implicit none

      integer, allocatable :: tnums(:)

      integer :: trtype

      logical :: skipID,cleantr,nobound

      ! number of tracks
      ntracks=0
      ncleantr=0
      ctrack=0

      write(*,*)"======================================="
      write(*,*)"=========== BUILDING TRACKS ==========="
      write(*,*)"======================================="

      ! open the file for tracks
      open(unit=1,file="tracks_all.txt",action="write",status="replace")
      open(unit=2,file="tracks_clean.txt",action="write",status="replace")

      ! now find connection trees for cells with nbw==0
      do i=1,globnIDs
        ! reset track type
        trtype=0
        cleantr=.false.
        nobound=.false.
        ! conditions for track initiation and track type
        skipID=.false.
        !if(tsclID(i)==1)skipID=.true.
        if(nbw(i)==0)then
          trtype=trtype+1
        else if(nbw(i)>1)then
          trtype=trtype+2
        else if(nbw(i)==1)then
          ! check if the bw link cell has more than one fw links -> cell splits from other track
          do k=1,nlinks(i)
            if(ltype(i,k)==1)then
              if(nfw(links(i,k))<2)then
                skipID=.true.
              else
                trtype=trtype+4             ! the code for this track initiation is 4
              end if
              exit
            end if
          end do
        else
          skipID=.true.
        end if

        if(skipID)cycle

        ntracks=ntracks+1
        ctrack=ctrack+1

        if(MOD(ctrack,outstep)==0 .OR. ctrack==1)then
          write(*,'(1a14,1i12)')" Found track: ",ctrack
        end if
        
        if(verbose)write(*,*)"Cell with ID ",clIDs(i)," initiates a track at ts: ",tsclID(i)
        
        ! allocate buffer for IDs and timesteps
        tp=1 ! this is the counter for the buffers
        allocate(tnums(maxtrlen))
        tnums=-1
        tnums(tp)=i

        k=i
        do
          ! the following ifs are the conditions for ending a track
          if(nfw(k)==0)then
            trtype=trtype+8
            ! backup if ...
            if(nbw(k)>1 .AND. tp>1)then
              tp=tp-1
              trtype=trtype-8+32
            end if
            exit
          else if(nbw(k)>1)then
            if(tp==1 .AND. trtype.ne.2)then
              trtype=trtype+32
              exit
            else if(tp==1 .AND. trtype==2 .AND. nfw(k)>1)then
              trtype=trtype+16
              exit
            else if(tp>1)then
              tp=tp-1
              trtype=trtype+32
              exit
            end if
          else if(nfw(k)>1)then
            trtype=trtype+16
            exit
          end if
          ! now walk through the tree of links
          do j=1,nlinks(k)
            if(ltype(k,j)==0)then ! take the first cell which is a forward link
              if(verbose)write(*,*)"continues at ts ",tsclID(j), " with ",clIDs(j)
              tp=tp+1
              tnums(tp)=links(k,j)
              k=links(k,j)
              exit
            end if
          end do
        end do

        ! is this a clean track?
        if(.NOT.ANY(touchb(tnums(1:tp))) .AND. tsclID(tnums(1)).ne.1 &
          &.AND. tsclID(tnums(tp)).ne.ntp)nobound=.true.

        if(trtype==9 .AND. nobound)then
          ncleantr=ncleantr+1
          cleantr=.true.
        end if

        ! if this track is clean/has no merges and splits write to tracks_clean.txt
        ! write header for this track
        if(cleantr)then
          write(2,'(1a4,1i12,1L4,1i4)')"### ",ctrack,nobound,trtype
          write(2,*)"       clID"!      tsclID      clarea    clcmassX    clcmassY"
          do k=1,tp
            write(2,'(1i12)')clIDs(tnums(k))!,tsclID(tnums(k)), &
              !& clarea(tnums(k)),clcmass(tnums(k),:)
          end do
        end if
        ! in any case write this track to file tracks_all.txt
        ! write header for this track
        write(1,'(1a4,1i12,1L4,1i4)')"### ",ctrack,nobound,trtype
        write(1,*)"       clID"!      tsclID      clarea    clcmassX    clcmassY"
        do k=1,tp
          write(1,'(1i12)')clIDs(tnums(k))!,tsclID(tnums(k)), &
            !& clarea(tnums(k)),clcmass(tnums(k),:)
        end do
        ! deallocate buffer
        deallocate(tnums)
      end do
      close(unit=1)
      close(unit=2)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,*)"Found ",ntracks," tracks ..."
      write(*,*)ncleantr," tracks do not split/merge or touch the boundaries."
      write(*,*)"---------"

      write(*,*)"======================================="
      write(*,*)"====== FINISHED BUILDING TRACKS ======="
      write(*,*)"======================================="

    end subroutine dobuildtracks

end module buildtracks
