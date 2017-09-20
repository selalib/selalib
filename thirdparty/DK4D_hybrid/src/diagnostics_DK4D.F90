!===========================================================================
!> Computation and printing of all the diagnostics
!>
!> \date 2014-08-25
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module diagnostics_DK4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_collective
  use equilibrium_DK4D_module
  use fdistribu_DK4D_module
  use field3d_DK4D_module
  use mesh_DK4D_module
  use sll_mpi, only : mpi_sum

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: diagnostics_DK4D_t

    !--> time associated to the diagnostics
    sll_real64 :: iter_time

    !--> Diagnostics for the norms
    sll_real64 :: mass
    sll_real64 :: L1_norm
    sll_real64 :: L2_norm
    sll_real64 :: Linf_norm
    sll_real64 :: entropy_kin

    !--> Diagnostics for energy
    sll_real64 :: energy_kin
    sll_real64 :: energy_pot
    sll_real64 :: energy_tot
    sll_real64 :: heat_flux

    !---> Diagnostiocs for Phi**2
    sll_real64 :: Phisquare

  end type diagnostics_DK4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Diagnostics: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_diagnostics_DK4D( diagnostics )

    type(diagnostics_DK4D_t), intent(inout) :: diagnostics

    sll_int32 :: ierr

  end subroutine new_diagnostics_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Diagnostics: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_diagnostics_DK4D( diagnostics )

    type(diagnostics_DK4D_t), intent(inout) :: diagnostics

    sll_int32 :: ierr

  end subroutine delete_diagnostics_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Diagnostics: Printing in HDF5 format
  !---------------------------------------------------------------------------
  subroutine print_diagnostics_DK4D( diagnostics, &
      idiag_num, &
      diag_time, &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      Phi )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
         sll_o_hdf5_ser_write_array, sll_s_hdf5_ser_file_close, &
         sll_t_hdf5_ser_handle

    type(diagnostics_DK4D_t), intent(inout) :: diagnostics
    sll_int32               , intent(in)    :: idiag_num
    sll_real64              , intent(in)    :: diag_time
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t), intent(in)    :: equilibrium
    type(fdistribu_DK4D_t)  , intent(in)    :: fdistribu
    type(field3d_DK4D_t)    , intent(in)    :: Phi

    !-> Local variables
    sll_int32 :: my_rank
    !---> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle
    character(len=80)    :: filename_HDF5
    character(20) , save :: numfmt = "'_d',i5.5"
    !---> For time
    sll_real64, dimension(1) :: diag_time_tmp
    sll_int32 , dimension(1) :: idiag_num_tmp
    !---> For norms
    sll_real64, dimension(1) :: mass_tmp
    sll_real64, dimension(1) :: L1_norm_tmp
    sll_real64, dimension(1) :: L2_norm_tmp
    sll_real64, dimension(1) :: Linf_norm_tmp
    sll_real64, dimension(1) :: entropy_kin_tmp
    !--> For energy
    sll_real64, dimension(1) :: energy_kin_tmp
    sll_real64, dimension(1) :: energy_pot_tmp
    !---> For Phi**2
    sll_real64, dimension(1) :: Phisquare_tmp

    write(filename_HDF5,'(A,'//numfmt//',A)') &
        'diags', idiag_num, ".h5"
    
    diag_time_tmp(1) = diag_time
    idiag_num_tmp(1) = idiag_num

    !*** Computation of the different diagnostics ***
    !--> Different norms
    call compute_L1_L2_Linf_norms( &
      mesh4d, &
      fdistribu, &
      diagnostics%mass, &
      diagnostics%L1_norm, &
      diagnostics%L2_norm, &
      diagnostics%Linf_norm, &
      diagnostics%entropy_kin )
    !--> Kinetic and potential energies
    call compute_energy( &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      Phi, &
      diagnostics%energy_kin, &
      diagnostics%energy_pot, &
      diagnostics%Phisquare )

    !*** Writing in HDF5 file ***
    !\todo write a routine in SELALIB to be able to save a real value
    mass_tmp(1)        = diagnostics%mass
    L1_norm_tmp(1)     = diagnostics%L1_norm
    L2_norm_tmp(1)     = diagnostics%L2_norm
    Linf_norm_tmp(1)   = diagnostics%Linf_norm
    entropy_kin_tmp(1) = diagnostics%entropy_kin
    energy_kin_tmp(1)  = diagnostics%energy_kin
    energy_pot_tmp(1)  = diagnostics%energy_pot
    Phisquare_tmp(1)   = diagnostics%Phisquare

    my_rank = sll_f_get_collective_rank(sll_v_world_collective)
    if (my_rank.eq.0) then
      call sll_s_hdf5_ser_file_create( trim(filename_HDF5), handle, file_err )
      call sll_o_hdf5_ser_write_array( handle, diag_time_tmp, &
          'diag_time', file_err )
      call sll_o_hdf5_ser_write_array( handle, mass_tmp, &
          'mass', file_err )
      call sll_o_hdf5_ser_write_array( handle, L1_norm_tmp, &
          'L1_norm', file_err )
      call sll_o_hdf5_ser_write_array( handle, L2_norm_tmp, &
          'L2_norm', file_err )
      call sll_o_hdf5_ser_write_array( handle, Linf_norm_tmp, &
          'Linf_norm', file_err )
      call sll_o_hdf5_ser_write_array( handle, entropy_kin_tmp, &
          'entropy_kin', file_err )
      call sll_o_hdf5_ser_write_array( handle, energy_kin_tmp, &
          'energy_kin', file_err )
      call sll_o_hdf5_ser_write_array( handle, energy_pot_tmp, &
          'energy_pot', file_err )
      call sll_o_hdf5_ser_write_array( handle, Phisquare_tmp, &
          'Phisquare', file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )

    end if

  end subroutine print_diagnostics_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Computation of the L1 norm , i.e   
  !   L1_norm = \int abs(delta f) * jac dvpar deta1 deta2 deta3
  !   delta f = f(x,y,z,vpar)
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !     
  !  Out: L1_norm
  !
  ! Computation of the L2 norm , i.e   
  !   L2_norm = sqrt(\int abs(delta f)**2 * jac dvpar deta1 deta2 deta3)
  !   delta f = f(x,y,z,vpar)
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !      
  !  Out: L2_norm
  !
  ! Computation of the L infini = max( abs( delta f) )
  !---------------------------------------------------------------------------
  subroutine compute_L1_L2_Linf_norms( &
      mesh4d, &
      fdistribu, &
      mass, &
      L1_norm, &
      L2_norm, &
      Linf_norm, &
      entropy_kin )

    use sll_m_collective

    type(mesh_DK4D_t)     , intent(in)  :: mesh4d
    type(fdistribu_DK4D_t), intent(in)  :: fdistribu
    sll_real64, intent(out) :: mass
    sll_real64, intent(out) :: L1_norm
    sll_real64, intent(out) :: L2_norm
    sll_real64, intent(out) :: Linf_norm
    sll_real64, intent(out) :: entropy_kin

    
    !-> Local variables
    !--> For local and global index 
    sll_int32  :: Neta1, Neta2, Neta3, Nvpar
    sll_int32  :: Neta1_loc, Neta2_loc
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1, i2, i3, i4
    sll_int32, dimension(1:4) :: glob_ind4d
    !--> For norm computations
    sll_real64 :: delta_eta1, delta_eta2, delta_eta3, delta_vpar
    sll_real64 :: jacob_tmp, val_jac, dphase_space
    sll_real64 :: f_tmp
    sll_real64 :: mass_loc
    sll_real64 :: L1_norm_loc, L2_norm_loc, Linf_norm_loc
    sll_real64 :: entropy_kin_loc
    !--> Temporary storage used for MPI
    sll_int32, parameter :: nb_diags = 5
    sll_real64, dimension(1:nb_diags) :: norm_diags_loc
    sll_real64, dimension(1:nb_diags) :: norm_diags

    Neta1_loc  = size(fdistribu%val4d_seqx3x4,1)
    Neta2_loc  = size(fdistribu%val4d_seqx3x4,2)
    Neta1      = mesh4d%eta1_eta2_mesh2d%num_cells1 + 1
    Neta2      = mesh4d%eta1_eta2_mesh2d%num_cells2 + 1
    Neta3      = size(fdistribu%val4d_seqx3x4,3)
    Nvpar      = size(fdistribu%val4d_seqx3x4,4)

    delta_eta1   = mesh4d%eta1_eta2_mesh2d%delta_eta1
    delta_eta2   = mesh4d%eta1_eta2_mesh2d%delta_eta2
    delta_eta3   = mesh4d%eta3_mesh1d%delta_eta
    delta_vpar   = mesh4d%vpar_mesh1d%delta_eta
    dphase_space = delta_eta1 * delta_eta2 * delta_eta3 * delta_vpar 

    L1_norm_loc     = 0.0_f64
    L2_norm_loc     = 0.0_f64
    Linf_norm_loc   = 0.0_f64
    mass_loc        = 0.0_f64
    entropy_kin_loc = 0.0_f64

    !*** Computation of the different diagnostics ***
    !***  locally in (x1,x2) directions           ***
    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        do i3 = 1,Neta3-1
          do i4 = 1,Nvpar-1

            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,i3,i4/) )           

            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)

            if ( i1 .ne. Neta1 ) then
              if ( i2 .ne. Neta2 ) then 

                jacob_tmp = mesh4d%transf_eta1eta2_xy%jacobian_at_node(i1,i2)
                val_jac   = abs(jacob_tmp)

                f_tmp = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4)

                mass_loc = mass_loc + f_tmp * val_jac * dphase_space

                L1_norm_loc = L1_norm_loc + &
                    abs(f_tmp) * val_jac * dphase_space

                L2_norm_loc = L2_norm_loc + &
                    sqrt( abs(f_tmp)**2 * val_jac * dphase_space )

                if ( f_tmp /= 0.0_f64 ) &
                    entropy_kin_loc = entropy_kin_loc - &
                    f_tmp * log(abs(f_tmp)) * val_jac * dphase_space

                Linf_norm_loc = max( abs(f_tmp), Linf_norm_loc )

              end if
            end if
          end do
        end do
      end do
    end do

    !*** Computation of the global diagnostics ***
    norm_diags_loc(1) = mass_loc
    norm_diags_loc(2) = L1_norm_loc
    norm_diags_loc(3) = L2_norm_loc
    norm_diags_loc(4) = Linf_norm_loc
    norm_diags_loc(5) = entropy_kin_loc

    call sll_s_collective_reduce_real64( &
        sll_v_world_collective, &
        norm_diags_loc, &
        nb_diags, &
        mpi_sum, &
        0, &
        norm_diags )

    mass        = norm_diags(1)
    L1_norm     = norm_diags(2)
    L2_norm     = norm_diags(3)
    Linf_norm   = norm_diags(4)
    entropy_kin = norm_diags(5)

  end subroutine compute_L1_L2_Linf_norms
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Computation of the kinetic energy, i.e   
  !   nrj_kin = \int delta f * vpar**2 * 0.5 * jac dvpar deta1 deta2 deta3
  !   delta f = f(x,y,z,vpar) - feq(x,y,vpar) 
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !  Out: nrj_kin
  !
  ! Computation of the potential energy, i.e   
  !   nrj_pot = \int delta f * Phi * 0.5 * jac dvpar deta1 deta2 deta3
  !   delta f = f(eta1,eta2,eta3,vpar) - feq(eta1,eta2,vpar) 
  !   Phi     = Phi(eta1,eta2,eta3) 
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !       Phi3d_seqx3(x1=distrib,x2=distrib,x3=*)
  !  Out: nrj_pot
  !
  ! Computation of the energy total = nrj_pot + nrj_kin
  !
  ! Computation of the heat flux, i.e   
  !   nrj_pot = \int delta f * Phi * 0.5 * jac dvpar deta1 deta2 deta3
  !   delta f = f(eta1,eta2,eta3,vpar) - feq(eta1,eta2,vpar) 
  !   Phi     = Phi(eta1,eta2,eta3) 
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !       Phi3d_seqx3(x1=distrib,x2=distrib,x3=*)
  !  Out: nrj_pot
  !
  !-----------------------------------------------------------  
  subroutine compute_energy( &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      Phi, &
      energy_kin, &
      energy_pot, &
      Phisquare )
      
    type(mesh_DK4D_t)       , intent(in) :: mesh4d
    type(equilibrium_DK4D_t), intent(in) :: equilibrium
    type(fdistribu_DK4D_t)  , intent(in) :: fdistribu
    type(field3d_DK4D_t)    , intent(in) :: Phi
    sll_real64, intent(out) :: energy_kin
    sll_real64, intent(out) :: energy_pot
    sll_real64, intent(out) :: Phisquare

    !-> Local variables
    !--> For local and global index 
    sll_int32  :: Neta1, Neta2, Neta3, Nvpar
    sll_int32  :: Neta1_loc, Neta2_loc
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1, i2, i3, i4
    sll_int32, dimension(1:4) :: glob_ind4d
    !--> For energy computations
    sll_real64 :: delta_eta1, delta_eta2, delta_eta3, delta_vpar
    sll_real64 :: jacob_tmp, val_jac
    sll_real64 :: dx_space, dphase_space
    sll_real64 :: vpar
    sll_real64 :: val_elec_pot
    sll_real64 :: delta_f
    sll_real64 :: energy_kin_loc, energy_pot_loc
    sll_real64 :: Phisquare_loc

    !--> Temporary storage used for MPI
    sll_int32, parameter :: nb_diags = 3
    sll_real64, dimension(1:nb_diags) :: energy_diags_loc
    sll_real64, dimension(1:nb_diags) :: energy_diags

    Neta1_loc  = size(fdistribu%val4d_seqx3x4,1)
    Neta2_loc  = size(fdistribu%val4d_seqx3x4,2)
    Neta1      = mesh4d%eta1_eta2_mesh2d%num_cells1 + 1
    Neta2      = mesh4d%eta1_eta2_mesh2d%num_cells2 + 1
    Neta3      = size(fdistribu%val4d_seqx3x4,3)
    Nvpar      = size(fdistribu%val4d_seqx3x4,4)

    delta_eta1   = mesh4d%eta1_eta2_mesh2d%delta_eta1
    delta_eta2   = mesh4d%eta1_eta2_mesh2d%delta_eta2
    delta_eta3   = mesh4d%eta3_mesh1d%delta_eta
    delta_vpar   = mesh4d%vpar_mesh1d%delta_eta
    dx_space     = delta_eta1 * delta_eta2 * delta_eta3
    dphase_space = dx_space * delta_vpar 

    energy_kin_loc = 0.0_f64
    energy_pot_loc = 0.0_f64
    Phisquare_loc  = 0.0_f64

    !*** Computation of the different diagnostics ***
    !***  locally in (x1,x2) directions           ***
    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        do i3 = 1,Neta3-1
          val_elec_pot = Phi%val3d_seqx3(iloc1,iloc2,i3)
          do i4 = 1,Nvpar-1
            vpar = mesh4d%vpar_grid(i4)

            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,i3,i4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)

            if (i1 .ne. Neta1) then
              if (i2 .ne. Neta2) then 

                jacob_tmp = mesh4d%transf_eta1eta2_xy%jacobian_at_node(i1,i2)
                val_jac   = abs(jacob_tmp)

                delta_f   = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4) - &
                    equilibrium%feq_xyvpar(i1,i2,i4)

                energy_kin_loc = energy_kin_loc + &
                    delta_f * vpar**2 * 0.5 * val_jac * dphase_space

                energy_pot_loc = energy_pot_loc + &
                    delta_f * val_elec_pot * 0.5 * val_jac * dphase_space

              end if
            end if
          end do

          if (i1 .ne. Neta1) then
            if (i2 .ne. Neta2) then 
              Phisquare_loc = Phisquare_loc + &
                  val_elec_pot**2 * val_jac * dphase_space 
            end if
          end if
        end do
      end do
    end do

    !*** Computation of the global diagnostics ***
    energy_diags_loc(1) = energy_kin_loc
    energy_diags_loc(2) = energy_pot_loc
    energy_diags_loc(3) = Phisquare_loc

    call sll_s_collective_reduce_real64( &
        sll_v_world_collective, &
        energy_diags_loc, &
        nb_diags, &
        mpi_sum, &
        0, &
        energy_diags )

    energy_kin = energy_diags(1)
    energy_pot = energy_diags(2)
    Phisquare  = energy_diags(3)

  end subroutine compute_energy
  !---------------------------------------------------------------------------

end module diagnostics_DK4D_module
!---------------------------------------------------------------------------
