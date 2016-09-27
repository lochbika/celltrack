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
module buildmetatracks

  use globvar

  implicit none

  contains
    subroutine dobuildmetatracks

      use globvar

      implicit none

      integer, allocatable :: trcon(:,:),conbuffer(:,:),metabuffer(:),metbuffer1(:),metbuffer2(:)
      integer :: maxtrcons,ncon
      real(kind=8) :: pint,avint
      logical :: nobound

      write(*,*)"======================================="
      write(*,*)"======== BUILDING META TRACKS ========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Searching for connections betw tracks ..."
      write(*,*)"---------"

      allocate(trcon(SUM(nfw),2))
      trcon=-1
      l=1
      do i=1,ntracks
        if(trtypes(i)==9)cycle
        ! get last cell of this track
        tp=0
        do k=1,maxtrlen
          if(alltracks(i,k)==-1)then
            tp=k-1
            exit
          else if(k==maxtrlen)then
            tp=k
            exit
          end if
        end do
        ! get its number in links matrix
        j=alltracks(i,tp)
        do k=i,ntracks
          tp=alltracks(k,1)
          if(abs(tsclID(tp)-tsclID(j))>1)cycle
          if(links(j,tp-minclIDloc(j)) .AND. i.ne.k)then
            trcon(l,1)=i
            trcon(l,2)=k
            if(MOD(l,outstep*2)==0 .OR. l==1)then
              write(*,*)"Found ",l," connections between tracks."
            end if
            l=l+1
          end if
        end do
      end do
      ncon=l-1
      write(*,*)"Found ",ncon," connections between tracks."
      ! this was the last time we needed links -> deallocate
      deallocate(links)

      write(*,*)"======================================="
      write(*,*)"=== Constructing meta tracks ..."
      write(*,*)"---------"

      open(unit=1,file="meta_all.txt",action="write",status="replace")
      open(unit=2,file="meta_con.txt",action="write",status="replace")

      ! what is the highest number of links between cells?
      if(MAXVAL(nfw)>MAXVAL(nbw))then
        maxtrcons=MAXVAL(nfw)*MAXVAL(nfw)*MAXVAL(nfw)
      else
        maxtrcons=MAXVAL(nbw)*MAXVAL(nbw)*MAXVAL(nbw)
      end if

      allocate(metabuffer(5000000),metbuffer1(maxtrcons),metbuffer2(maxtrcons))
      allocate(conbuffer(SUM(nfw),2))
      cmeta=0
      nmeta=0

      do
        if(ALL(trcon==-1))exit
        metabuffer=-1
        conbuffer=-1
        metbuffer1=-1
        tp=1
        l=1
        do i=1,ncon
          if(ANY(trcon(i,:)==-1))cycle
          ! initiate a meta track by types 17 and 33
          if(trtypes(trcon(i,1))==17 .OR. trtypes(trcon(i,1))==33)then
            ! initiate the buffers and counters
            metabuffer(tp)=trcon(i,1)
            metabuffer(tp+1)=trcon(i,2)
            metbuffer1(1)=trcon(i,1)
            metbuffer1(2)=trcon(i,2)
            conbuffer(1,:)=trcon(i,:)
            l=l+1
            trcon(i,:)=-1
            tp=tp+2
            cmeta=cmeta+1
            if(MOD(cmeta,outstep/2)==0 .OR. cmeta==1)then
              write(*,'(1a19,1i6)')" Found meta track: ",cmeta
            end if
            exit
          end if
        end do

        do
          if(ALL(metbuffer1==-1))exit
          k=1
          metbuffer2=-1
          do i=1,maxtrcons
            if(metbuffer1(i)==-1)exit
            if(.NOT.ANY(metabuffer==metbuffer1(i)))then
              metabuffer(tp)=metbuffer1(i)
              tp=tp+1
            end if
          end do
          do j=1,maxtrcons
            if(metbuffer1(j)==-1)exit
            do i=1,ncon
              if(trcon(i,1)==metbuffer1(j))then
                metbuffer2(k)=trcon(i,2)
                conbuffer(l,:)=trcon(i,:)
                trcon(i,:)=-1
                k=k+1
                l=l+1
              end if
              if(trcon(i,2)==metbuffer1(j))then
                metbuffer2(k)=trcon(i,1)
                conbuffer(l,:)=trcon(i,:)
                trcon(i,:)=-1
                k=k+1
                l=l+1
              end if
            end do
          end do
          metbuffer1=metbuffer2
        end do
        ! save this meta track to file
        ! check if any of the tracks touches the boundaries
        nobound=.true.
        do k=1,maxtrcons
          if(metabuffer(k)==-1)exit
          if(.NOT.nobounds(metabuffer(k)))nobound=.false.
        end do

        ! write header meta_con.txt
        write(2,'(1a4,1i12,1L4,1i4)')"### ",cmeta,nobound
        write(2,*)"   trackID1    trackID2"
        do i=1,ncon
          if(conbuffer(i,1)==-1)exit
          write(2,'(2i12)')conbuffer(i,1), conbuffer(i,2)
        end do

        ! write header meta_all.txt
        write(1,'(1a4,1i12,1L4,1i4)')"### ",cmeta,nobound
        write(1,*)"    trackID"

        ! write tack details
        do k=1,maxtrcons
          if(metabuffer(k)==-1)exit
          ! calc peak value
          pint=0
          do j=1,maxtrlen
            if(alltracks(metabuffer(k),j)==-1)exit
            if(pint<clpint(alltracks(metabuffer(k),j)))pint=clpint(alltracks(metabuffer(k),j))
          end do
          ! calc average value
          avint=0
          tp=0
          do j=1,maxtrlen
            if(alltracks(metabuffer(k),j)==-1)exit
            avint=avint+clavint(alltracks(metabuffer(k),j))
            tp=tp+1
          end do
          avint=avint/tp
          write(1,'(1i12)')metabuffer(k)!,trtypes(metabuffer(k)),pint,avint
        end do

      end do
      nmeta=cmeta

      close(unit=1)
      close(unit=2)

      deallocate(trcon,conbuffer,metbuffer1,metbuffer2,metabuffer)

      write(*,*)"======================================="
      write(*,*)"=== Summary ..."
      write(*,*)"Found ",nmeta," meta tracks."
      write(*,*)"---------"

      write(*,*)"======================================="
      write(*,*)"==== FINISHED BUILDING META TRACKS ===="
      write(*,*)"======================================="

    end subroutine dobuildmetatracks

end module buildmetatracks
