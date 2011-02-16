!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!------------------------------------------------------
! file : read_profile.f90
! date : 10/09/2009
!  Reading of radial profiles given by input files
!-----------------------------------------------------
module read_profile_module
  use prec_const
  use globals, only : pglobal_id
  use geometry_class
      
  implicit none
  include "mpiperso.h"
      
  !******************************
  contains
  !******************************
      
  !---------------------------------------------------- 
  !---------------------------------------------------- 
  subroutine read_ascii_profile(profile_filename,geom,interp_profile)
    use globals     , only : rhomin, rhomax, a
    use utils_module, only : locate
    use spline1d_class
    character(LEN=*) , intent(in) :: profile_filename  
    type(geometry)   , intent(in) :: geom
    real(RKIND), &
        dimension(:) , pointer    :: interp_profile   ! (0:geom%Nr)
      
    integer     :: uprof = 30
    integer     :: ios, ieof, ierr
    integer     :: icount, ipoint, Nbpoints, ir, nbelem
    real(RKIND) :: x_tmp, y_tmp
    !-> for data profile reading
    real(RKIND)                        :: dr_prof
    real(RKIND), dimension(:), pointer :: data_pos, data_prof
    !-> for cubic spline interpolation
    !----> Hermite boundary conditions 
    integer          , parameter :: BCr_left  = 2  
    integer          , parameter :: BCr_right = 2
    type(nspline1d)              :: nspline1d_r
    integer                      :: i, ipos
    real(RKIND), dimension(-1:2) :: sbase_prof
    real(RKIND)                  :: rho_interp, fct_interp
      
    if (pglobal_id.eq.0) then
      open(uprof,file=trim(profile_filename), &
        status="old", iostat=ios)
      !*** count the number of lines of the file ***
      icount = 0
      if (ios.eq.0) then
        read(uprof,*,iostat=ieof) x_tmp, y_tmp
        do while (ieof==0)
          icount = icount + 1
          read(uprof,*,iostat=ieof) x_tmp, y_tmp
        end do
      else 
        print*,'The file ',profile_filename, ' does not exist'
        stop
      end if
      Nbpoints = icount-1
      if (Nbpoints.lt.0) then
        print*, 'Problem in the reading of the file ', &
          profile_filename,' => Nbpoints = ',Nbpoints
        stop
      else
        print*,'Number of points read in ', &
          profile_filename,' = ',Nbpoints
      end if
      
      !*** read the radial position of the points + ***
      !***   the profile associated                 ***
      rewind(unit=uprof)      
      allocate(data_pos(0:Nbpoints))
      allocate(data_prof(0:Nbpoints))
      do ipoint = 0,Nbpoints
        read(uprof,*) data_pos(ipoint), data_prof(ipoint)
      end do
      close(uprof)
      
      !*** test if the positions of the points are normalized, ***
      !***  i.e between 0 and 1                                ***
      if ( (data_pos(0).lt.0._RKIND) .or. &
        (data_pos(0).gt.rhomin) ) then
        print*,' The first point must be between 0 and rhomin = ', &
          rhomin
        stop
      end if
      
      if ( (data_pos(Nbpoints).lt.rhomax) .or. &
        (data_pos(Nbpoints).gt.1._RKIND) ) then
        print*,' The last point must be between rhomax = ', &
          rhomax,' and 1'
        stop
      end if
      
      !*** interpolation to compute the profile ***
      !***  on the radial grid                  ***
      dr_prof = abs(data_pos(1)-data_pos(0))
      !-> computation of the cubic spline coefficients     
      call new_spline1d_natural(nspline1d_r,Nbpoints,dr_prof)
      call natural_spline_coef(nspline1d_r,data_prof, &
        BCr_left,BCr_right)
      !-> interpolation to compute the value of the profile 
      !->   on the radial mesh grid
      do ir = 0,geom%Nr
        rho_interp = geom%rg(ir)/a
        call locate(rho_interp,data_pos,Nbpoints,dr_prof, &
          " read_ascii_profile "//char(0),ipos)
        call spline_basis(data_pos(ipos),rho_interp, &
          data_pos(ipos+1),dr_prof,sbase_prof)
        fct_interp = 0._RKIND
        do i = -1,2
          fct_interp = fct_interp + &
            nspline1d_r%scoef(ipos+i)*sbase_prof(i)
        end do
        interp_profile(ir) = fct_interp
      end do
      
      call del_spline1d_natural(nspline1d_r)
      deallocate(data_pos)
      deallocate(data_prof)
    end if
      
    !*** distribution of the profile on all the processors ***
    nbelem = geom%Nr + 1
      
    CALL MPI_BCAST(interp_profile,nbelem,MPI_REAL8, &
      0,MPI_COMM_WORLD,ierr)
  end subroutine read_ascii_profile
end module read_profile_module
