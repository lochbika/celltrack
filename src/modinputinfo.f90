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
module inputinfo

  implicit none

  contains
  
    subroutine gatherInputInfo()
      
      use ncdfpars
    
      use globvar, only : ifile,ntp,status,tp,ivar,stdclen,vunit,vname,vdate,vtime, &
        & tstep,nx,ny,minx,maxx,miny,maxy,level,inmissval,nlev,xunit,yunit,diflon,diflat,thresdir,thres
    
      write(*,*)"======================================="
      write(*,*)"========== INPUT INFORMATION =========="
      write(*,*)"======================================="
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to input file..."
    
      ! Get Information about the variable
      CALL getVarInfo(ifile)
    
      ! Get initial Information about time axis
      CALL getTaxisInfo(ifile)

      ! Get grid information
      CALL gethorGrid(ifile)
      
      ! Get Zaxis info
      CALL getZaxisInfo(ifile)
      
      ! print input summary to stdout
      write(*,*)"======================================="
      write(*,*)"=== INPUT SUMMARY:"
      write(*,*)"---------"
      write(*,*)"Input   : ",trim(ifile)
      write(*,*)"---------"
      write(*,'(A,1a12)')" VAR        : ",trim(vname)
      write(*,'(A,1a12)')" Unit       : ",trim(vunit)
      if(inmissval.eq.-123456789.D0)then
        write(*,'(A,1f12.2)')" MissVal    :   undefined!"
      else
        write(*,'(A,1f12.2)')" MissVal    : ",inmissval
      end if
      write(*,'(A,1i12)')" NX         : ",nx
      write(*,'(A,1f12.2)')" MIN X      : ",minx
      write(*,'(A,1f12.2)')" MAX X      : ",maxx
      write(*,'(A,1f12.2)')" DIF X      : ",diflon
      write(*,'(A,1a12)')" Unit       : ",trim(xunit)
      write(*,'(A,1i12)')" NY         : ",ny
      write(*,'(A,1f12.2)')" MIN Y      : ",miny
      write(*,'(A,1f12.2)')" MAX Y      : ",maxy
      write(*,'(A,1f12.2)')" DIF Y      : ",diflat
      write(*,'(A,1a12)')" Unit       : ",trim(yunit)
      write(*,'(A,1i8,1i8.6)')" START DATE : ",vdate(1),vtime(1)
      write(*,'(A,1i8,1I8.6)')" END DATE   : ",vdate(ntp),vtime(ntp)
      write(*,'(A,1i12)')" NTSTEPS    : ",ntp
      write(*,'(A,1i12)')" TSTEP      : ",tstep
      write(*,'(A,1i12)')" NLEV       : ",nlev
      write(*,'(A,1f12.2)')" SELLEV     : ",level
      write(*,*)"---------"
      if(thresdir.eq.1)then
        write(*,'(A,1f12.2)')" threshold : MIN ",thres
      else
        write(*,'(A,1f12.2)')" threshold : MAX ",thres
      end if
      write(*,*)"======================================="
      write(*,*)"======== END INPUT INFORMATION ========"
      write(*,*)"======================================="
    
    end subroutine gatherInputInfo
    
end module inputinfo
