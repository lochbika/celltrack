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
module mainstreamdetection

  use globvar

  implicit none

  contains
    subroutine domainstreamdetection

      use globvar

      implicit none

      integer :: ant,run,maxconlen,numeqsol
      integer, allocatable :: init(:),quit(:),qi(:),ninit,nquit
      integer, allocatable :: paths(:,:),qicon(:,:)
      integer, allocatable :: thismainstream(:),lastmainstream(:)
      real(kind=8),allocatable :: eta(:),pher(:)
      logical,allocatable :: backw(:)
      logical :: resetnants
      character(len=800) :: ttrack
      real(kind=8) :: wsum,isum,areasum,intsum

      write(*,*)"======================================="
      write(*,*)"======== MAINSTREAM DETECTION ========="
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Preparing ..."
      write(*,*)"---------"

      ! allocate mainstream
      allocate(allmainstream(nmeta,maxmetalen))
      allmainstream=-1

      ! read connections file
      CALL scan_trackfile("meta_con.txt",nmeta,maxconlen)
      allocate(allcon(nmeta,maxconlen,2))
      allcon=-1

      ! open the connections file
      open(unit=1,file="meta_con.txt",action="read",status="old")

      !now loop all lines
      cmeta=0

      do
        ! read the line ;)
        read(1,'(A)',IOSTAT=riostat)ttrack
        if(riostat.ne.0)exit
        !check if a new track begins
        if(ttrack(1:3)=="###")then
          read(1,*)
          ! buffer counter
          tp=1
          read(ttrack(5:16),'(1i12)')cmeta
          if(verbose)then
            if(MOD(cmeta,outstep)==0 .OR. cmeta==1 .OR. cmeta==nmeta)then
              write(*,*)"Reading connections for meta track ",cmeta," of ",nmeta
            end if
          end if
        else
          !now read cell ID, time step and area
          read(ttrack,'(2i12)')allcon(cmeta,tp,:)
          tp=tp+1
        end if
      end do

      close(unit=1)

      write(*,*)"======================================="
      write(*,*)"=== Do it! ..."
      write(*,*)"---------"

      ! open file for connections and pheromone values
      open(unit=1,file="meta_con_pher.txt",action="write",status="replace")

      ! loop meta tracks
      do i=1,nmeta
        ! status to stdout
        if(MOD(i,outstep/2)==0 .OR. i==1 .OR. i==nmeta)then
          write(*,*)"Processing meta track ",i," of ",nmeta
        end if

        ! allocate and init array and counter for stagnation detection
        allocate(thismainstream(maxmetalen),lastmainstream(maxmetalen))
        thismainstream=-1
        lastmainstream=-1
        numeqsol=0

        ! which tracks are of type 17 and 33 -> init
        allocate(init(maxmetalen))
        init=-1
        tp=0
        do k=1,maxmetalen
          if(allmeta(i,k)==-1)exit
          if(trtypes(allmeta(i,k))==17 .OR. trtypes(allmeta(i,k))==33)then
            tp=tp+1
            init(tp)=allmeta(i,k)
          end if
        end do
        ninit=tp
        ! which tracks are of type 10 and 12 -> quit
        allocate(quit(maxmetalen))
        quit=-1
        tp=0
        do k=1,maxmetalen
          if(allmeta(i,k)==-1)exit
          if(trtypes(allmeta(i,k))==10 .OR. trtypes(allmeta(i,k))==12)then
            tp=tp+1
            quit(tp)=allmeta(i,k)
          end if
        end do
        nquit=tp
        ! if not set, set the number of ants equal to the number of connections
        if(nants.eq.-1)then
          resetnants=.true.
          do k=1,maxconlen
            if(allcon(i,k,1)==-1)exit
            nants=k
          end do
          if(nants<5)nants=5
          if(nants>maxnants .AND. maxnants.ne.-1)nants=maxnants
          if(verbose)write(*,*)"--- number of ants: ",nants
        else
          resetnants=.false.
        end if
        ! calculate the step of peakVAL and clarea for all connections and initialize eta
        allocate(eta(maxconlen))
        eta=-1
        do k=1,maxconlen
          if(allcon(i,k,1)==-1)exit
          ! get last cell of trackID1
          do j=1,maxtrlen
            if(alltracks(allcon(i,k,1),j)==-1)exit
            intsum=clpint(alltracks(allcon(i,k,1),j))
            areasum=clarea(alltracks(allcon(i,k,1),j))
          end do
          intsum=sqrt( ( intsum - clpint(alltracks(allcon(i,k,2),1)) )**2 ) /intsum
          areasum=sqrt( ( areasum - clarea(alltracks(allcon(i,k,2),1)) )**2 ) / areasum
          eta(k)=(intsum + areasum )/2
          ! backup; if the distance between cells is 0
          if(eta(k)==0.D0)then
            if(verbose)write(*,*)"Warning: difference between cells is 0!!"
            eta(k)=0.000001
          end if
        end do
        ! initialize the pheromone tracks
        allocate(pher(maxconlen))
        pher=-1
        allocate(paths(nants,maxconlen*2))
        ! do this for each ant to determine more trustworthy initial pheromone values
        isum=0 ! total average cost value
        do n=1,nants
          paths=-1
          ! use a random init and do a nn track; temporarely use paths to store it
          CALL RANDOM_NUMBER(rnum)
          tp=1 + FLOOR( rnum * ninit )
          ! construct a nn path
          CALL nnPath(init(tp),allcon(i,:,:),maxconlen,eta,paths(1,:))
          ! calculate the average cost (sum of eta) of this random path
          wsum=0
          tp=0
          do k=1,(maxconlen*2) !path
            if(paths(1,k)==-1)exit
            do j=1,maxconlen !connections
              if(allcon(i,j,1)==-1)exit
              if(allcon(i,j,1)==paths(1,k) .AND. allcon(i,j,2)==paths(1,k+1))then
                wsum=wsum+eta(j)
                tp=tp+1
              end if
            end do
          end do
          ! average
          wsum=wsum/tp
          ! sum up total wsum
          isum=isum+wsum
        end do
        ! average isum and save in wsum
        wsum=isum/nants
        ! reset paths
        paths=-1
        ! update the pheromone trails with
        do k=1,maxconlen !pher
          if(allcon(i,k,1)==-1)exit
          pher(k)=nants / wsum
        end do

        ! now release the ants ;)
        do run=1,nruns

          ! allocate the backwards marker
          allocate(backw(nants))
          backw=.false.

          do ant=1,nants

            ! get a random init position
            CALL RANDOM_NUMBER(rnum)
            tp=1 + FLOOR( rnum * ( ninit+nquit) )
            ! check if tp is either in the init or quit range and set qi to init resp. quit
            allocate(qicon(maxconlen,2))
            if(tp>ninit)then
              allocate(qi(nquit))
              qi=quit(1:nquit)
              tp=tp-ninit
              backw(ant)=.true.
              ! swap the connection direction
              do k=1,maxconlen
                if(allcon(i,k,1)==-1)exit
                qicon(k,1)=allcon(i,k,2)
                qicon(k,2)=allcon(i,k,1)
              end do
            else
              allocate(qi(ninit))
              qi=init(1:ninit)
              qicon=allcon(i,:,:)
              backw(ant)=.false.
            end if

            ! construct a path
            CALL acoPath(qi(tp),qicon,maxconlen,eta,pher,paths(ant,:))

            ! deallocate qicon,qi
            deallocate(qicon,qi)
          end do

          ! pheromone evaporation
          do k=1,maxconlen
            if(pher(k)==-1)exit
            pher(k) = (1-pherevap)*pher(k)
          end do

          ! pheromone update
          do ant=1,nants

            ! calculate the cost (sum of eta) of each ants path
            wsum=0
            tp=0
            do k=1,(maxconlen*2) !path
              if(paths(ant,k)==-1)exit
              do j=1,maxconlen !connections
                if(allcon(i,j,1)==-1)exit
                if(backw(ant))then
                  if(allcon(i,j,2)==paths(ant,k) .AND. allcon(i,j,1)==paths(ant,k+1))then
                    wsum=wsum+eta(j)
                    tp=tp+1
                  end if
                else
                  if(allcon(i,j,1)==paths(ant,k) .AND. allcon(i,j,2)==paths(ant,k+1))then
                    wsum=wsum+eta(j)
                    tp=tp+1
                  end if
                end if
              end do
            end do

            ! average by number of nodes
            wsum=wsum/tp

            ! update the pheromone path
            do k=1,maxconlen ! connections
              if(allcon(i,k,1)==-1)exit
              ! search a path for this connection
              do j=1,(maxconlen*2)
                if(paths(ant,j)==-1)exit
                if(backw(ant))then
                  if(allcon(i,k,2)==paths(ant,j) .AND. allcon(i,k,1)==paths(ant,j+1))then
                    pher(k)=pher(k)+1/wsum
                  end if
                else
                  if(allcon(i,k,1)==paths(ant,j) .AND. allcon(i,k,2)==paths(ant,j+1))then
                    pher(k)=pher(k)+1/wsum
                  end if
                end if
              end do
            end do
          end do
          deallocate(backw)

          ! Do the following each run
          if(MOD(run,1)==0 .OR. run==nruns)then
            ! check if the current mainstream is different from the last ones
            ! generate the solution: chose the init with the highest pheromone trail and do nn path with the inverse pheromone values
            wsum=0
            tp=0
            do k=1,maxconlen
              if(allcon(i,k,1)==-1)exit
              if(ANY(allcon(i,k,1)==init))then
                if(pher(k)>wsum)then
                  wsum=pher(k)
                  tp=k
                end if
              end if
            end do
            ! nn path
            thismainstream=-1
            CALL nnPath(allcon(i,tp,1),allcon(i,:,:),maxconlen,1/pher,thismainstream(:))
            if(ALL(thismainstream==lastmainstream))then
              numeqsol=numeqsol+1
            else
              numeqsol=0
            end if
            lastmainstream=thismainstream
          end if

          ! check the termination conditions
          if(numeqsol>10)then
            if(verbose)write(*,*)"Stagnation in solution construction: Mainstream didn't change since 10 iterations. Stop!"
            exit
          end if
          if(run==nruns)then
            ! only notify if in verbose mode; the loop will stop after this anyway
            if(verbose)write(*,*)"Maximum number of runs reached. Stop!"
          end if
        end do

        ! use the latest thismainstream to set the global mainstream for this meta track
        allmainstream(i,:)=thismainstream

        ! write the connections with their pheromone values and distances
        ! write header meta_con_pher.txt
        write(1,'(1a4,1i12,1L4,1i4)')"### ",i,mnobounds(i)
        write(1,*)"   trackID1    trackID2        pher        dist"
        do k=1,maxmetalen
          if(allcon(i,k,1)==-1)exit
          write(1,'(2i12,2f12.5)')allcon(i,k,:),pher(k),eta(k)
        end do

        deallocate(init,quit,eta,pher,paths,thismainstream,lastmainstream)
        if(resetnants)nants=-1
      end do
      close(unit=1)

      open(unit=1,file="meta_mainstream.txt",action="write",status="replace")
      do i=1,nmeta
        ! write header meta_stats.txt
        write(1,'(1a4,1i12,1L4,1i4)')"### ",i,mnobounds(i)
        write(1,*)"    trackID      trType            peakVal    pValtime              avVal start   dur"
        do k=1,maxmetalen
          if(allmainstream(i,k)==-1)exit
          write(1,'(2i12,1f19.12,1i12,1f19.12,2i6)')allmainstream(i,k),trtypes(allmainstream(i,k)), &
          & trpint(allmainstream(i,k)),trpinttime(allmainstream(i,k)), &
          & travint(allmainstream(i,k)),tsclID(alltracks(allmainstream(i,k),1)),trdur(allmainstream(i,k))
        end do
      end do

      ! save to which meta track mainstream a cell belongs
      allocate(clmetamstr(globnIDs))
      clmetamstr=-1
      do i=1,nmeta
        do k=1,maxmetalen
          if(allmainstream(i,k)==-1)exit
          do j=1,maxtrlen
            if(alltracks(allmainstream(i,k),j)==-1)exit
            clmetamstr(alltracks(allmainstream(i,k),j))=i
          end do
        end do
      end do

      write(*,*)"======================================="
      write(*,*)"==== FINISHED MAINSTREAM DETECTION ===="
      write(*,*)"======================================="

    end subroutine domainstreamdetection

    subroutine acoPath(init,cons,ncons,lens,pher,path)

      use globvar, only: alpha,beta

      implicit none

      integer,intent(out) :: path(ncons*2)
      integer,intent(in)  :: cons(ncons,2),ncons,init
      integer :: a,next(500),pcount,i,tp
      real(kind=8),intent(in) :: pher(ncons),lens(ncons)
      real(kind=8) :: tauij,etaij,tauil,etail,zeta,rnum,probsum
      real(kind=8),allocatable :: cprob(:)

      pcount=0
      a=init
      do
        ! add a to path
        pcount=pcount+1
        path(pcount)=a
        ! search for possible connections at this cell
        tp=0 ! number of possible connections
        next=-1
        do i=1,ncons
          if(cons(i,1)==-1)exit
          if(cons(i,1)==a)then
            tp=tp+1
            next(tp)=i
          end if
        end do
        ! there are no possible connections? exit this loop
        if(ALL(next==-1))exit
        !!! the random proportional rule; but do it only if there are at least 2 possible choices
        if(tp==1)then
          a=cons(next(tp),2)
        else
          ! calculate zeta
          zeta=0
          do i=1,tp
            tauil=pher(next(i))
            etail=1/lens(next(i))
            zeta=zeta+( tauil**alpha * etail**beta )
          end do
          ! calc tauij and etaij and the probability of decisions
          allocate(cprob(tp))
          do i=1,tp
            tauij=pher(next(i))
            etaij=1/lens(next(i))
            cprob(i)=( tauij**alpha * etaij**beta ) / zeta
          end do
          ! generate a random number
          CALL RANDOM_NUMBER(rnum)
          probsum=0
          do i=1,tp
            probsum=probsum+cprob(i)
            if(rnum<=probsum)then
              a=cons(next(i),2)
              exit
            end if
          end do
          deallocate(cprob)
        end if
      end do
    end subroutine acoPath

    subroutine nnPath(init,cons,ncons,lens,path)
      !! if lens contains differences this algorithm will find a nearest neighbour path
      !! if lens contains absolute values this is a local minimum cost algorithm

      implicit none

      integer,intent(out) :: path(ncons*2)
      integer,intent(in)  :: cons(ncons,2),ncons,init

      real(kind=8),intent(in) :: lens(ncons)

      real(kind=8) :: clen

      integer :: i,tp,a,next(500),pcount

      pcount=0
      a=init
      do
        ! add a to path
        pcount=pcount+1
        path(pcount)=a
        ! search for possible connections at this cell
        tp=0
        next=-1
        do i=1,ncons
          if(cons(i,1)==-1)exit
          if(cons(i,1)==a)then
            tp=tp+1
            next(tp)=i
          end if
        end do
        ! there are no possible connections? exit this loop
        if(ALL(next==-1))exit
        ! now choose the nearest neighbor
        clen=HUGE(clen)
        do i=1,500
          if(next(i)==-1)exit
          if(lens(next(i))<clen)then
            clen=lens(next(i))
            a=cons(next(i),2)
          end if
        end do
      end do

    end subroutine nnPath

end module mainstreamdetection
