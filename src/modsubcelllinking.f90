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
module subcelllinking

  use globvar

  implicit none


  contains
    subroutine dosubcelllinking

      use globvar
      use ncdfpars

      implicit none

      include 'cdi.inc'

      ! data arrays
      real(kind=stdfloattype), allocatable :: subcells(:)     ! array for reading from subcells nc
      real(kind=stdfloattype), allocatable :: cells(:)        ! array for reading from cells nc

      write(*,*)"======================================="
      write(*,*)"======== START SUBCELL LINKING ========"
      write(*,*)"======================================="

      write(*,*)"======================================="
      write(*,*)"=== Opening connection to cells file..."

      ! Open the dataset
      streamID1=streamOpenRead(outfile)
      if(streamID1<0)then
         write(*,*)cdiStringError(streamID1)
         stop
      end if

      ! Set the variable IDs 1
      varID1=0
      vlistID1=streamInqVlist(streamID1)
      gridID1=vlistInqVarGrid(vlistID1,varID1)
      taxisID1=vlistInqTaxis(vlistID1)
      zaxisID1=vlistInqVarZaxis(vlistID1,varID1)
      
      write(*,*)"======================================="
      write(*,*)"=== Opening connection to SUBcells file..."

      ! Open the dataset
      streamID2=streamOpenRead(suboutfile)
      if(streamID2<0)then
         write(*,*)cdiStringError(streamID2)
         stop
      end if

      ! Set the variable IDs 1
      varID2=0
      vlistID2=streamInqVlist(streamID2)
      gridID2=vlistInqVarGrid(vlistID2,varID2)
      taxisID2=vlistInqTaxis(vlistID2)
      zaxisID2=vlistInqVarZaxis(vlistID2,varID2)
      
      ! allocate and initialize linking array
      allocate(sublinks(globsubnIDs))
      sublinks=-1

      ! Get data
      write(*,*)"======================================="
      write(*,*)"=== Linking SUBcells to cells ..."
      do tsID=0,(ntp-1)
        if(MOD(tsID+1,outstep)==0 .OR. tsID==0 .OR. tsID==ntp-1 .OR. verbose)then
          write(*,*)"Processing timestep: ",tsID+1,"/",ntp,"..."
        end if

        ! Allocate arrays for data storage
        allocate(subcells(nx*ny))
        allocate(cells(nx*ny))

        ! Set time step for input file
        status=streamInqTimestep(streamID1,tsID)

        ! Set time step for cells file
        status=streamInqTimestep(streamID2,tsID)

        ! Read time step from cells input
        call streamReadVarSlice(streamID1,varID1,0,cells,nmiss1)

        ! Read time step from SUBcells input
        call streamReadVarSlice(streamID2,varID2,0,subcells,nmiss2)
        
        ! cycle if field contains only missing values; but write it to output
        if(nmiss1==nx*ny)then
          deallocate(subcells,cells)
          cycle
        end if

        ! no iterate through the input array and link subcells to cells
        do i=1,nx*ny
          if(cells(i)==outmissval)then
            cycle
          else
            if(sublinks(INT(subcells(i)))==-1)sublinks(INT(subcells(i)))=INT(cells(i))
          end if
        end do
        
        ! clean up
        deallocate(subcells,cells)
      end do
      
      ! close input and output
      CALL streamClose(streamID1)
      CALL streamClose(streamID2)

      write(*,*)"======================================="
      write(*,*)"====== FINISHED SUBCELL LINKING ======="
      write(*,*)"======================================="

    end subroutine dosubcelllinking

end module subcelllinking
