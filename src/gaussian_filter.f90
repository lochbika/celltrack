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
!-----
!! originally from https://github.com/nichannah/gaussian-filter/blob/master/gaussian_filter.F90

subroutine gaussian_kernel(sigma, truncate)

    use globvar, only: kernel,stdfloattype

    real(kind=stdfloattype), intent(in) :: sigma
    integer, intent(in), optional :: truncate

    real(kind=stdfloattype), allocatable :: x(:,:), y(:,:)
    integer :: radius, trunc, i, j
    real(kind=stdfloattype) :: s

    if (present(truncate)) then
        trunc = truncate
    else
        trunc = 4
    endif

    radius = NINT(trunc * sigma + 0.5)
    s = sigma**2

    ! Set up meshgrid.
    allocate(x(-radius:radius, -radius:radius))
    allocate(y(-radius:radius, -radius:radius))
    do j = -radius, radius
        do i = -radius, radius
            x(i, j) = i
            y(i, j) = j
        enddo
    enddo

    ! Make kernel; deallocate first if done before
    if(allocated(kernel)) then
      deallocate(kernel)
    endif
    
    allocate(kernel(-radius:radius, -radius:radius))
    kernel = 2.0*exp(-0.5 * (x**2 + y**2) / s)
    kernel = kernel / sum(kernel)

    deallocate(x)
    deallocate(y)

end subroutine gaussian_kernel
