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
      
!-----------------------------------------------------------------
! file : coord_system.f90
! date : 07/17/2009
!  Initialisation of the system of coordinates (x1,x2,x3),
!   i.e :
!    - metric tensor defined as the following matrix:
!        g11 g12 g13
!        g21 g22 g23
!        g31 g32 g33
!    - jacobian in space
!    - norm of a vector in the system
!-----------------------------------------------------------------
module coord_system_class
  use prec_const
  use globals, only : Nr, Ntheta
  use mem_alloc_module
  use geometry_class
      
  implicit none
      
  !*** coordinate_system definition ***
  !-> 1/R(r,theta) definition 
  !     with R(r,theta) = RO+h(r)*cos(theta)  
  real(RKIND), dimension(:,:), pointer, public :: R
  real(RKIND), dimension(:,:), pointer, public :: inv_R
      
  type :: coord_system
    !-> components of the metric tensor
    real(RKIND), dimension(:,:), pointer :: g11
    real(RKIND), dimension(:,:), pointer :: g12
    real(RKIND), dimension(:,:), pointer :: g13
    real(RKIND), dimension(:,:), pointer :: g21
    real(RKIND), dimension(:,:), pointer :: g22
    real(RKIND), dimension(:,:), pointer :: g23
    real(RKIND), dimension(:,:), pointer :: g31
    real(RKIND), dimension(:,:), pointer :: g32
    real(RKIND), dimension(:,:), pointer :: g33
    !> components of the contravariant metric tensor
    real(RKIND), dimension(:,:), pointer :: gradx1gradx1
    real(RKIND), dimension(:,:), pointer :: gradx1gradx2
    real(RKIND), dimension(:,:), pointer :: gradx1gradx3
    real(RKIND), dimension(:,:), pointer :: gradx2gradx1
    real(RKIND), dimension(:,:), pointer :: gradx2gradx2
    real(RKIND), dimension(:,:), pointer :: gradx2gradx3
    real(RKIND), dimension(:,:), pointer :: gradx3gradx1
    real(RKIND), dimension(:,:), pointer :: gradx3gradx2
    real(RKIND), dimension(:,:), pointer :: gradx3gradx3
  end type coord_system
 
  !-> jacobian in the real space
  real(RKIND), dimension(:,:), pointer, public :: jacobian_space
  real(RKIND), dimension(:)  , pointer, public :: intdtheta_Js
  real(RKIND), dimension(:)  , pointer, public :: intdthetadphi_Js
      
  !******************************
  contains
  !******************************
     
  !---------------------------------------------------- 
  ! Constructor 
  !----------------------------------------------------
  subroutine new_coord_system(cthis)
    use globals, only : Nr, Ntheta
    type(coord_system), intent(out) :: cthis
      
    !-> allocation for R(r,theta)
    call glob_allocate(R,0,Nr,0,Ntheta,'R')
    call glob_allocate(inv_R,0,Nr,0,Ntheta,'inv_R')
    
    !-> allocation for the metric tensor covariant components
    call glob_allocate(cthis%g11,0,Nr,0,Ntheta,'g11')
    call glob_allocate(cthis%g12,0,Nr,0,Ntheta,'g12')
    call glob_allocate(cthis%g13,0,Nr,0,Ntheta,'g13')
    call glob_allocate(cthis%g21,0,Nr,0,Ntheta,'g21')
    call glob_allocate(cthis%g22,0,Nr,0,Ntheta,'g22')
    call glob_allocate(cthis%g23,0,Nr,0,Ntheta,'g23')
    call glob_allocate(cthis%g31,0,Nr,0,Ntheta,'g31')
    call glob_allocate(cthis%g32,0,Nr,0,Ntheta,'g32')
    call glob_allocate(cthis%g33,0,Nr,0,Ntheta,'g33')
      
    !-> allocation for the metric tensor contravariant components
    call glob_allocate(cthis%gradx1gradx1, &
      0,Nr,0,Ntheta,'gradx1gradx1')
    call glob_allocate(cthis%gradx1gradx2, &
      0,Nr,0,Ntheta,'gradx1gradx2')
    call glob_allocate(cthis%gradx1gradx3, &
      0,Nr,0,Ntheta,'gradx1gradx3')
    call glob_allocate(cthis%gradx2gradx1, &
      0,Nr,0,Ntheta,'gradx2gradx1')
    call glob_allocate(cthis%gradx2gradx2, &
      0,Nr,0,Ntheta,'gradx2gradx2')
    call glob_allocate(cthis%gradx2gradx3, &
      0,Nr,0,Ntheta,'gradx2gradx3')
    call glob_allocate(cthis%gradx3gradx1, &
      0,Nr,0,Ntheta,'gradx3gradx1')
    call glob_allocate(cthis%gradx3gradx2, &
      0,Nr,0,Ntheta,'gradx3gradx2')
    call glob_allocate(cthis%gradx3gradx3, &
      0,Nr,0,Ntheta,'gradx3gradx3')
      
    !-> allocation for the jacobian in space
    call glob_allocate(jacobian_space, &
      0,Nr,0,Ntheta,'jacobian_space')
    !-> \int dtheta Js with Js the jacobian in space
    call glob_allocate(intdtheta_Js,0,Nr,'intdtheta_Js')
    !-> \int dtheta dphi Js with Js the jacobian in space
    call glob_allocate(intdthetadphi_Js,0,Nr,'intdthetadphi_Js')
  end subroutine new_coord_system
  
      
  !---------------------------------------------------- 
  ! destructor
  !---------------------------------------------------- 
  subroutine del_coord_system(cthis)
    type(coord_system), intent(inout) :: cthis
      
    !-> deallocation for R(r,theta)
    call glob_deallocate(R)
    call glob_deallocate(inv_R)
      
    !-> deallocation for the metric tensor covariant components
    call glob_deallocate(cthis%g11)
    call glob_deallocate(cthis%g12)
    call glob_deallocate(cthis%g13)
    call glob_deallocate(cthis%g21)
    call glob_deallocate(cthis%g22)
    call glob_deallocate(cthis%g23)
    call glob_deallocate(cthis%g31)
    call glob_deallocate(cthis%g32)
    call glob_deallocate(cthis%g33)
      
    !-> deallocation for the metric tensor contravariant components
    call glob_deallocate(cthis%gradx1gradx1)
    call glob_deallocate(cthis%gradx1gradx2)
    call glob_deallocate(cthis%gradx1gradx3)
    call glob_deallocate(cthis%gradx2gradx1)
    call glob_deallocate(cthis%gradx2gradx2)
    call glob_deallocate(cthis%gradx2gradx3)
    call glob_deallocate(cthis%gradx3gradx1)
    call glob_deallocate(cthis%gradx3gradx2)
    call glob_deallocate(cthis%gradx3gradx3)
      
    !-> deallocation for the jacobian in space
    call glob_deallocate(jacobian_space)
    call glob_deallocate(intdtheta_Js)
    call glob_deallocate(intdthetadphi_Js)
  end subroutine del_coord_system
      
  !---------------------------------------------------- 
  ! - initialization of  R(r,theta) = R0+h(r)*cos(theta)
  !  where h(r) = r for the moment
  !---------------------------------------------------- 
  subroutine init_R(geom)
    use globals, only : slab_geometry, R0
    type(geometry), intent(in) :: geom
    
    real(RKIND) :: hri
    integer     :: ir, itheta
      
    !*** Initialization of 1/R ***
    if (.not.slab_geometry) then
      do ir = 0,Nr
        hri = geom%rg(ir)
        do itheta = 0,Ntheta
          !-> R = R0+h(r)*cos(theta)
          R(ir,itheta) = R0 + hri*geom%cos_theta(itheta)
          !-> 1/R = 1/(R0+h(r)*cos(theta))
          inv_R(ir,itheta) = 1._RKIND/R(ir,itheta)
        end do
      end do
    else
      !*** R=R0 in the gradients for the 4D case ***
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          R(ir,itheta)     = R0
          inv_R(ir,itheta) = 1._RKIND/R(ir,itheta)
        end do
      end do
    end if
  end subroutine init_R
      
  !---------------------------------------------------- 
  ! Initialisation of the covariant components
  !  of the metric tensor
  !----------------------------------------------------
  subroutine init_metric_tensor(cthis,geom)
    type(coord_system) , intent(out) :: cthis
    type(geometry)     , intent(in)  :: geom
      
    integer     :: ir, itheta
    real(RKIND) :: ri, Rij
      
    do itheta = 0,Ntheta
      do ir = 0,Nr
        ri  = geom%rg(ir)
        Rij = R(ir,itheta)
        cthis%g11(ir,itheta) = 1._RKIND
        cthis%g12(ir,itheta) = 0._RKIND
        cthis%g13(ir,itheta) = 0._RKIND
        cthis%g21(ir,itheta) = 0._RKIND
        cthis%g22(ir,itheta) = ri*ri
        cthis%g23(ir,itheta) = 0._RKIND
        cthis%g31(ir,itheta) = 0._RKIND
        cthis%g32(ir,itheta) = 0._RKIND
        cthis%g33(ir,itheta) = Rij*Rij       
      end do
    end do
  end subroutine init_metric_tensor
      
  !---------------------------------------------------- 
  ! Computation of the jacobian
  !  Jacobian_space = sqrt(g) with
  !   g = sqrt{g11*g22*g33-g11*g32*g23-g12*g21*g33
  !            +g12*g31*g23+g13*g21*g32-g13*g31*g22}
  !----------------------------------------------------
  subroutine compute_jacobian_space(cthis,geom)
    use globals, only : Nr, Ntheta, Nphi
    type(coord_system), intent(in) :: cthis
    type(geometry)    , intent(in) :: geom
      
    integer     :: ir, itheta, iphi
    real(RKIND) :: g11, g12, g13, g21, g22
    real(RKIND) :: g23, g31, g32, g33
    real(RKIND) :: jacob_space_tmp    
      
    do itheta = 0,Ntheta
      do ir = 0,Nr
        g11 = cthis%g11(ir,itheta)
        g12 = cthis%g12(ir,itheta)
        g13 = cthis%g13(ir,itheta)
        g21 = cthis%g21(ir,itheta)
        g22 = cthis%g22(ir,itheta)
        g23 = cthis%g23(ir,itheta)
        g31 = cthis%g31(ir,itheta)
        g32 = cthis%g32(ir,itheta)
        g33 = cthis%g33(ir,itheta)     
        jacobian_space(ir,itheta) = &
          sqrt(g11*g22*g33 - g11*g32*g23 - &
          g12*g21*g33 + g12*g31*g23 + &
          g13*g21*g32 - g13*g31*g22)
      end do
    end do
    !-> \int dtheta Js and \int dtheta dphi Js 
    !->    with Js is the jacobian in space
    do ir = 0,Nr
      intdtheta_Js(ir)     = 0._RKIND
      intdthetadphi_Js(ir) = 0._RKIND
      do itheta = 0,Ntheta
        jacob_space_tmp  = jacobian_space(ir,itheta)
        intdtheta_Js(ir) = intdtheta_Js(ir) + &
          jacob_space_tmp*geom%coeff_intdtheta(itheta)
        do iphi = 0,Nphi
          intdthetadphi_Js(ir) = intdthetadphi_Js(ir) + &
            jacob_space_tmp*geom%coeff_intdtheta(itheta)* &
            geom%coeff_intdphi(iphi)        
        end do
      end do
    end do
  end subroutine compute_jacobian_space
      
  !---------------------------------------------------- 
  ! Computation of the elements of the contravariant
  !  metric tensor {g^ij}=grad_xi.grad_xj as:
  !    g22g33-g23g32  g13g32-g12g33  g12g23-g22g13
  !    g23g31-g21g33  g11g33-g13g31  g13g21-g11g23
  !    g21g32-g22g31  g12g31-g11g32  g11g22-g21g12
  !----------------------------------------------------
  subroutine compute_contravariant_mtensor(cthis)
    use globals, only : Nr, Ntheta
    type(coord_system), intent(inout) :: cthis
      
    integer     :: ir, itheta
    real(RKIND) :: g11, g12, g13, g21, g22
    real(RKIND) :: g23, g31, g32, g33    
    real(RKIND) :: jacob_tmp, inv_jacob2
    
    do itheta = 0,Ntheta
      do ir = 0,Nr
        jacob_tmp  = jacobian_space(ir,itheta)
        inv_jacob2 = 1._RKIND/(jacob_tmp*jacob_tmp)
        g11        = cthis%g11(ir,itheta)
        g12        = cthis%g12(ir,itheta)
        g13        = cthis%g13(ir,itheta)
        g21        = cthis%g21(ir,itheta)
        g22        = cthis%g22(ir,itheta)
        g23        = cthis%g23(ir,itheta)
        g31        = cthis%g31(ir,itheta)
        g32        = cthis%g32(ir,itheta)
        g33        = cthis%g33(ir,itheta)     
        cthis%gradx1gradx1(ir,itheta) = (g22*g33-g23*g32) * &
          inv_jacob2
        cthis%gradx1gradx2(ir,itheta) = (g13*g32-g12*g33) * &
          inv_jacob2
        cthis%gradx1gradx3(ir,itheta) = (g12*g23-g22*g13) * &
          inv_jacob2
        cthis%gradx2gradx1(ir,itheta) = (g23*g31-g21*g33) * &
          inv_jacob2
        cthis%gradx2gradx2(ir,itheta) = (g11*g33-g13*g31) * &
          inv_jacob2
        cthis%gradx2gradx3(ir,itheta) = (g13*g21-g11*g23) * &
          inv_jacob2
        cthis%gradx3gradx1(ir,itheta) = (g21*g32-g22*g31) * &
          inv_jacob2
        cthis%gradx3gradx2(ir,itheta) = (g12*g31-g11*g32) * &
          inv_jacob2
        cthis%gradx3gradx3(ir,itheta) = (g11*g22-g21*g12) * &
          inv_jacob2
      end do
    end do
  end subroutine compute_contravariant_mtensor
      
  !------------------------------------------------------------
  !  Initialisation of the coordinate system
  !------------------------------------------------------------
  subroutine init_coordinate_system(cthis,geom)
    use globals, only : memory_test
    type(coord_system), intent(inout) :: cthis
    type(geometry)    , intent(in)    :: geom    
      
    !*** arrays allocations ***
    call new_coord_system(cthis)
      
    !*** initialisation of the coordinate system ***
    if (.not. memory_test) then
      !*** initialisation of R(r,theta) ***
      call init_R(geom)
      !*** initialisation of the metric tensor ***
      call init_metric_tensor(cthis,geom)
      !*** computation the jacobian in space ***
      call compute_jacobian_space(cthis,geom)
      !*** computation of the contravariant metric tensor ***
      call compute_contravariant_mtensor(cthis)
    end if
  end subroutine init_coordinate_system
      
  !------------------------------------------------------------
  ! Computation of the contravariant components of 
  !  a vector A by using its covariant components, i.e:
  !   -> input  : A1, A2 and A3
  !   -> output :
  !       . A_gradx1 = A1 gradx1.gradx1 + 
  !                    A2 gradx2.gradx1 + A3 gradx3.gradx1
  !       . A_gradx2 = A1 gradx1.gradx2 + 
  !                    A2 gradx2.gradx2 + A3 gradx3.gradx2
  !       . A_gradx3 = A1 gradx1.gradx3 + 
  !                    A2 gradx2.gradx3 + A3 gradx3.gradx3
  !------------------------------------------------------------
  subroutine compute_contravariant_vector(cthis,A1,A2,A3, &
    A_gradx1,A_gradx2,A_gradx3)
    type(coord_system)           , intent(in)  :: cthis
    real(RKIND), dimension(0:,0:), intent(in)  :: A1, A2, A3
    real(RKIND), dimension(0:,0:), intent(out) :: A_gradx1
    real(RKIND), dimension(0:,0:), intent(out) :: A_gradx2
    real(RKIND), dimension(0:,0:), intent(out) :: A_gradx3
      
    integer :: ir, itheta
    real(RKIND) :: A1_ij, A2_ij, A3_ij
    real(RKIND) :: gradx1gradx1_ij, gradx1gradx2_ij
    real(RKIND) :: gradx1gradx3_ij, gradx2gradx1_ij
    real(RKIND) :: gradx2gradx2_ij, gradx2gradx3_ij
    real(RKIND) :: gradx3gradx1_ij, gradx3gradx2_ij
    real(RKIND) :: gradx3gradx3_ij
      
    do itheta = 0,Ntheta
      do ir = 0,Nr
        A1_ij           = A1(ir,itheta)
        A2_ij           = A2(ir,itheta)
        A3_ij           = A3(ir,itheta)
        gradx1gradx1_ij = cthis%gradx1gradx1(ir,itheta)
        gradx1gradx2_ij = cthis%gradx1gradx2(ir,itheta)
        gradx1gradx3_ij = cthis%gradx1gradx3(ir,itheta)
        gradx2gradx1_ij = cthis%gradx2gradx1(ir,itheta)
        gradx2gradx2_ij = cthis%gradx2gradx2(ir,itheta)
        gradx2gradx3_ij = cthis%gradx2gradx3(ir,itheta)
        gradx3gradx1_ij = cthis%gradx3gradx1(ir,itheta)
        gradx3gradx2_ij = cthis%gradx3gradx2(ir,itheta)
        gradx3gradx3_ij = cthis%gradx3gradx3(ir,itheta)
        A_gradx1(ir,itheta) = A1_ij*gradx1gradx1_ij + & 
          A2_ij*gradx2gradx1_ij + A3_ij*gradx3gradx1_ij
        A_gradx2(ir,itheta) = A1_ij*gradx1gradx2_ij + &
          A2_ij*gradx2gradx2_ij + A3_ij*gradx3gradx2_ij
        A_gradx3(ir,itheta) = A1_ij*gradx1gradx3_ij + &
          A2_ij*gradx2gradx3_ij + A3_ij*gradx3gradx3_ij 
      end do
    end do
  end subroutine compute_contravariant_vector
      
  !---------------------------------------------------- 
  ! Computation of the norm of a vector A:
  !  -> input : covariant components of vec_A, i.e :
  !         A1, A2 and A3
  !  -> output :  A_norm = \sqrt{A1^2*g11
  !                        + A^2*g22 + A_3^2*g33}
  !----------------------------------------------------
  subroutine compute_norm(cthis,A1,A2,A3,Anorm)
    type(coord_system)           , intent(in)  :: cthis
    real(RKIND), dimension(0:,0:), intent(in)  :: A1
    real(RKIND), dimension(0:,0:), intent(in)  :: A2
    real(RKIND), dimension(0:,0:), intent(in)  :: A3
    real(RKIND), dimension(0:,0:), intent(out) :: Anorm
      
    integer     :: ir, itheta
    real(RKIND) :: A1_ij, A2_ij, A3_ij
    real(RKIND) :: gradx1gradx1_ij
    real(RKIND) :: gradx2gradx2_ij
    real(RKIND) :: gradx3gradx3_ij
      
    do itheta = 0,Ntheta
      do ir = 0,Nr
        A1_ij           = A1(ir,itheta)
        A2_ij           = A2(ir,itheta)
        A3_ij           = A3(ir,itheta)
        gradx1gradx1_ij = cthis%gradx1gradx1(ir,itheta)
        gradx2gradx2_ij = cthis%gradx2gradx2(ir,itheta)
        gradx3gradx3_ij = cthis%gradx3gradx3(ir,itheta)
        Anorm(ir,itheta) =  sqrt( &
          A1_ij*A1_ij*gradx1gradx1_ij + &
          A2_ij*A2_ij*gradx2gradx2_ij + &
          A3_ij*A3_ij*gradx3gradx3_ij )
      end do
    end do
  end subroutine compute_norm
end module coord_system_class
