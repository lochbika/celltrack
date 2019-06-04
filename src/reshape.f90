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
subroutine reshapeF2d(datin,nx,ny,datout)

  use globvar, only: stdfloattype

  implicit none
  integer                   :: tp,nx,ny,x,y
  real(kind=stdfloattype), intent(in)  :: datin(nx,ny)
  real(kind=stdfloattype), intent(out) :: datout(nx*ny)

  tp=1
  do y=1,ny
    do x=1,nx
      datout(tp)=datin(x,y)
      tp=tp+1
    end do
  end do

end subroutine reshapeF2d

subroutine reshapeT2d(datin,nx,ny,datout)

  use globvar, only: stdfloattype

  implicit none
  integer                   :: tp,nx,ny,x,y
  real(kind=stdfloattype), intent(in)  :: datin(nx*ny)
  real(kind=stdfloattype), intent(out) :: datout(nx,ny)

  tp=1
  do y=1,ny
    do x=1,nx
      datout(x,y)=datin(tp)
      tp=tp+1
    end do
  end do

end subroutine reshapeT2d
