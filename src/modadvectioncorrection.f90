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
module advectioncorrection

  use globvar

  implicit none

  contains
    subroutine doadvectioncorrection

      ! modules
      use globvar
      use ncdfpars
      use celldetection
      use cellstatistics
      use celllinking
      use linkstatistics

      implicit none

      ! data arrays
      real(kind=8), allocatable :: uvfield2d(:,:),uvfield(:),vvfield2d(:,:),vvfield(:)   ! the velocity fields in 1d and 2d
      integer :: vnx,vny,coarsex=12,coarsey=12,selCL
      real(kind=8), allocatable :: vxvals(:),vyvals(:),vclx(:),vcly(:)
      integer, allocatable :: vclxindex(:),vclyindex(:)      ! the indices of the neares gridpoints for each cell (x and y coord)
      integer, allocatable :: smpsize(:,:)                   ! sample size for each gridpoint on the velocity field
      real(kind=8) :: mindist,cdist

      write(*,*)"======================================="
      write(*,*)"===== START ADVECTION CORRECTION ======"
      write(*,*)"======================================="
      
      do adviter=1,nadviter
      
        write(*,*)"======================================="
        write(*,*)"=== This is iteration ",adviter
  
        ! at first do the cell detection and statistics
        CALL docelldetection()
        CALL calccellstatistics()
  
        ! now we do the linking and the statistics
        CALL linking()
        CALL calclinkstatistics()
  
        ! now we can use the gathered information to coarse grain the grid and open a new dataset
        vnx=nx/coarsex
        vny=ny/coarsey
        allocate(vxvals(vnx),vyvals(vny))
        do x=1,vnx
          vxvals(x)=(xvals(0)-diflon/2) + (diflon*coarsex*x) - (diflon*coarsex/2)
          write(*,*)vxvals(x)
        end do
        do y=1,vny
          vyvals(y)=(yvals(0)-diflat/2) + (diflat*coarsey*y) - (diflat*coarsey/2)
          write(*,*)vyvals(y)
        end do
  
        ! find the nearest gridpoint on the velocity grid for all cells
        allocate(vclxindex(globnIDs),vclyindex(globnIDs))
        do clID=1,globnIDs
          mindist=HUGE(mindist)
          do x=1,vnx
            cdist=abs(vxvals(x)-wclcmass(clID,1)*diflon)
            if(cdist<mindist)then
              vclxindex(clID)=x
              mindist=cdist
            end if
          end do
          mindist=HUGE(mindist)
          do y=1,vny
            cdist=abs(vyvals(y)-wclcmass(clID,2)*diflat)
            if(cdist<mindist)then
              vclyindex(clID)=y
              mindist=cdist
            end if
          end do
          !write(*,*)vclxindex(clID),vclyindex(clID)
        end do
  
        ! now calculate each cells velocity
        allocate(vclx(globnIDs),vcly(globnIDs))
        vclx=-999.D0
        vcly=-999.D0
        do clID=1,globnIDs
          if(tsclID(clID).ne.1 .AND. .NOT.touchb(clID))then
            if(nbw(clID)==1)then
              ! find the cell which is connected backwards
              do i=1,iclIDloc(clID)
                if(links(clID,i))then
                  selCL=i
                  exit
                end if
              end do
              vclx(clID)=(wclcmass(clID,1)-wclcmass(selCL,1))/30
              vcly(clID)=(wclcmass(clID,2)-wclcmass(selCL,2))/30
            end if
          end if
          !write(*,*)vclx(clID),vcly(clID)
        end do
  
        ! calculate average velocities on the grid for each time step
        do tsID=0,(ntp-1)
          allocate(uvfield2d(vnx,vny),vvfield2d(vnx,vny),smpsize(vnx,vny))
          uvfield2d=0
          vvfield2d=0
          smpsize=0
          
          do clID=1,globnIDs
            if(tsclID(clID)>tsID+1)exit
            if(tsclID(clID)==tsID+1)then
              uvfield2d(vclxindex(clID),vclyindex(clID)) = uvfield2d(vclxindex(clID),vclyindex(clID)) + vclx(clID)
              vvfield2d(vclxindex(clID),vclyindex(clID)) = vvfield2d(vclxindex(clID),vclyindex(clID)) + vcly(clID)
              smpsize(vclxindex(clID),vclyindex(clID)) = smpsize(vclxindex(clID),vclyindex(clID)) + 1
            end if
          end do
          
          ! average and set 0 sized gridpoints to missing value
          WHERE(smpsize.ne.0)uvfield2d=uvfield2d/smpsize
          WHERE(smpsize.ne.0)vvfield2d=vvfield2d/smpsize
          WHERE(smpsize==0)uvfield2d=-999.D0
          WHERE(smpsize==0)vvfield2d=-999.D0
          
          deallocate(uvfield2d,vvfield2d,smpsize)
        end do
        
        ! deallocate all arrays to rerun the detection and statistics part
        deallocate(xvals,yvals,levels,clids,tsclID,clarea,touchb,clcmass,wclcmass,clavint,clpint)
        deallocate(links,minclIDloc,iclIDloc,nbw,nfw)
        deallocate(vxvals,vyvals,vclxindex,vclyindex,vclx,vcly)
        
      end do
      
      write(*,*)"======================================="
      write(*,*)"==== FINISHED ADVECTION CORRECTION ===="
      write(*,*)"======================================="

    end subroutine doadvectioncorrection

end module advectioncorrection
