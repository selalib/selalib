!File: Module Poisson
!
! Poisson solver on unstructured mesh with triangles
!
module sll_m_poisson_2d_tri
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_choleski, only: &
      sll_s_choles, &
      sll_s_desrem

   use sll_m_triangular_meshes, only: &
      sll_t_triangular_mesh_2d

   implicit none

   public :: &
      sll_s_poisson_2d_triangular_init, &
      sll_s_poisson_2d_triangular_free, &
      sll_s_compute_e_from_phi, &
      sll_s_compute_e_from_rho, &
      sll_s_compute_phi_from_rho, &
      sll_t_poisson_2d_triangular

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   interface sll_o_create
      module procedure initialize_poisson_solver
      module procedure initialize_poisson_solver_from_file
   end interface sll_o_create

!    Mesh parameters
!
!        nbs     - number of nodes
!        nbt     - number of triangles
!        ndiric  - number of nodes Dirichlet
!        ndirb3  - number of nodes Dirichlet for Ampere
!        nbfrax  - number of nodes axis/boundary
!        nmxfr   - number of boundary referenced
!        nmxsd   - number of subdomains referenced
!        nefro   - number of elements with one node on the boundary
!        nelfr   - number of elements on the boundaries
!        nelin   - number of elements not on the boundaries
!        irefdir - boundary references Dirichlet non homogeneous
!        nnoeuf  - number of nodes Dirichlet non homogeneous
!
!    Geometry :
!
!        xlml   - limit below x of domain
!        xlmu   - limit above x of domain
!        ylml   - limit below y of domain
!        ylmu   - limit above y of domain
!        petitl - small length reference
!        grandl - big length reference
!        imxref - undefined reference
!
!    Delaunay-Voronoi:
!
!        nbcoti - number of edges Delaunay inside
!        nbtcot - number of edges Delaunay total
!        nbcfli - number of edges inside with wrong CFL
!        ncotcu - number of edges fo each boundary
!        nnref  - number of references Dirichlet non homogeneous boundary
!
!    Edges on boundaries:
!
!        nctfrt - number total de edges boundaries
!        nctfro - number de edges boundaries par reference
!        nctfrp - pointer des arrays de edges boundaries
!
!    Internal boundaries:
!
!        nnofnt - number nodes internal boundaries
!        ntrfnt - number total triangles
!        ntrfrn - number triangles par boundaries
!        ntrfrc - number cumsum triangles
!        nndfnt - nodes Dirichlet internal boundaries
!        ncdfnt - edges Dirichlet internal boundaries
!
!    Array elements types:
!
!        nmxcol - number of type
!        nclcol - number of elements in every type
!
   sll_real64, parameter :: grandx = 1.d+20

!>@brief
!> Derived type for Poisson solver on unstructured mesh with triangles.
!> @details
!> We using here P1-conformin finite element method.
!> This program is derived from F. Assous and J. Segre code M2V.
   type :: sll_t_poisson_2d_triangular

      private
      sll_real64, dimension(:), allocatable :: vnx    !< normal vector on node
      sll_real64, dimension(:), allocatable :: vny    !< normal vector on node
      logical, dimension(:), allocatable :: naux   !< true if the node is on boundary
      sll_real64, dimension(:), allocatable :: gradx  !< x-grad matrix
      sll_real64, dimension(:), allocatable :: grady  !< y-grad matrix
      sll_real64, dimension(:), allocatable :: grgr   !< grad-grad matrix
      sll_int32, dimension(:), allocatable :: mors1  !< array to store matrix
      sll_int32, dimension(:), allocatable :: mors2  !< array to store matrix
      sll_int32, dimension(:), allocatable :: iprof  !< array to store matrix
      sll_int32, dimension(:), allocatable :: ifron
      sll_real64, dimension(:), allocatable :: amass
      sll_real64, dimension(:), allocatable :: vtantx
      sll_real64, dimension(:), allocatable :: vtanty
      sll_real64, dimension(:), allocatable :: sv1
      sll_real64, dimension(:), allocatable :: sv2

      sll_int32  :: ndiric
      sll_int32  :: ntypfr(5)
      sll_real64 :: potfr(5)
      sll_real64 :: eps0 = 1.0_f64 !8.8542e-12

      type(sll_t_triangular_mesh_2d), pointer :: mesh => null()

   end type sll_t_poisson_2d_triangular

contains

   subroutine sll_s_poisson_2d_triangular_init(solver, mesh, ntypfr, potfr)

      type(sll_t_poisson_2d_triangular)             :: solver
      type(sll_t_triangular_mesh_2d), intent(in)  :: mesh
      sll_int32, dimension(:), intent(in)  :: ntypfr
      sll_real64, dimension(:), intent(in)  :: potfr

      call initialize_poisson_solver(solver, mesh, ntypfr, potfr)

   end subroutine sll_s_poisson_2d_triangular_init

!> Delete the solver derived type
   subroutine sll_s_poisson_2d_triangular_free(self)
      type(sll_t_poisson_2d_triangular) :: self

#ifdef DEBUG
      print *, 'delete poisson solver'
#endif

   end subroutine sll_s_poisson_2d_triangular_free

!> Compute electric field from charge density
!> @param[in]  self solver derived type
!> @param[in]  rho  charge density array on nodes
!> @param[out] phi  electric potential on nodes
!> @param[out] ex   electric field x component on nodes
!> @param[out] ey   electric field y component on nodes
   subroutine sll_s_compute_e_from_rho(self, rho, phi, ex, ey)
      type(sll_t_poisson_2d_triangular) :: self
      sll_real64, intent(in)          :: rho(:)
      sll_real64, intent(out)         :: phi(:)
      sll_real64, intent(out)         :: ex(:)
      sll_real64, intent(out)         :: ey(:)

      call poissn(self, rho, phi, ex, ey)
      call poifrc(self, ex, ey)

   end subroutine sll_s_compute_e_from_rho

   subroutine sll_s_compute_phi_from_rho(self, rho, phi)
      type(sll_t_poisson_2d_triangular) :: self
      sll_real64, intent(in)          :: rho(:)
      sll_real64, intent(out)         :: phi(:)

      call poissn(self, rho, phi)

   end subroutine sll_s_compute_phi_from_rho

   subroutine sll_s_compute_e_from_phi(self, phi, ex, ey)
      type(sll_t_poisson_2d_triangular) :: self
      sll_real64, intent(in)          :: phi(:)
      sll_real64, intent(out)         :: ex(:)
      sll_real64, intent(out)         :: ey(:)

      call poliss(self, phi, ex, ey)
      call poifrc(self, ex, ey)

   end subroutine sll_s_compute_e_from_phi

   subroutine read_data_solver(ntypfr, potfr)

      sll_int32, parameter       :: nfrmx = 5
      sll_int32, intent(out)     :: ntypfr(nfrmx)
      sll_real64, intent(out)     :: potfr(nfrmx)

      sll_int32                   :: ityp
      sll_int32                   :: ifr
      sll_int32                   :: i
      character(len=72)           :: argv
      character(len=132)          :: inpfil
      logical                     :: lask
      character(len=*), parameter :: self_sub_name = 'read_data_solver'

      NAMELIST /nlcham/ ntypfr, potfr

      call get_command_argument(1, argv); write (*, '(1x, a)') argv

!------------------------------------------------------------!
!     Reads in default parameters from input file (.inp)        !
!     If LASK=t, ask user if file is to be read.                !
!------------------------------------------------------------!

      lask = .true.

      inpfil = trim(argv)//'.inp'

      write (*, "(/10x,'Fichier d''entree : ',a)") inpfil

      open (10, file=inpfil, status='OLD', err=80)
      write (*, 1050, advance='no') trim(inpfil)
      lask = .false.

80    continue

      if (lask) then
         write (*, 1900) inpfil
         write (*, 1800, advance='no')
         read (*, "(a)") inpfil
         write (*, 1700) trim(inpfil)
      end if

      write (6, 900)

      ntypfr = 1
      potfr = 0.0_f64

      open (10, file=inpfil)
      write (*, *) "Lecture de la namelist nlcham"
      read (10, nlcham, err=100)
      goto 200
100   continue
      write (*, *) "Error de lecture de la namelist nlcham"
      write (*, *) "Valeurs par defaut"
      write (*, nlcham)
      stop
200   close (10)

      write (6, 921)

      do i = 1, nfrmx
         write (6, 922) i, ntypfr(i), potfr(i)
      end do

      write (6, 932)
      do ifr = 1, nfrmx
         ityp = ntypfr(ifr)
         if (ityp == 1) then
            write (6, 902) ifr
         else if (ityp == 3) then
            write (6, 904) ifr
         else
            write (6, 910) ifr, ityp
            SLL_ERROR(self_sub_name, "Unspecified error.")
         end if
      end do

900   format(//10x, 'Boundary conditions')
902   format(/10x, 'Reference :', i3, 5x, 'Dirichlet ')
904   format(/10x, 'Reference :', i3, 5x, 'Neumann')
910   format(/10x, 'Option not available'/  &
              &       10x, 'Reference :', i3, 5x, 'Type :', i3/)
921   format(//5x, 'Frontiere', 5x, 'ntypfr', 7x, 'potfr')
922   format(5x, I3, 4x, I10, 2E12.3, 2I10, E12.3)
932   format(//10x, 'Boundary conditions'/)
1050  format(/' Read settings from file  ', A, ' ?  Y')
1700  format(/' New parameters write to file  ', A,/)
1800  format(/' Settings may have been changed - New title :')
1900  format(/' Input file  ', A, '  not found')

   end subroutine read_data_solver

   subroutine initialize_poisson_solver_from_file(self, mesh)

      type(sll_t_poisson_2d_triangular), intent(out)     :: self
      type(sll_t_triangular_mesh_2d), intent(in), target :: mesh
      sll_int32                                         :: ntypfr(5)
      sll_real64                                        :: potfr(5)

      call read_data_solver(ntypfr, potfr)

      call initialize_poisson_solver(self, mesh, ntypfr, potfr)

   end subroutine initialize_poisson_solver_from_file

!Subroutine: init_solveur_poisson
! Allocate arrays
!
! amass - matrix mass diagonal
! mors1 - matrix morse2
! mors2 - elements matrices morse
! prof  - matrix profil
! grgr  - matrix grad-grad
! gradx - matrix gradx
! grady - matrix grady
! fron  - nodes Dirichlet
! xaux  - arrays auxiliar reals
! iaux  - arrays auxiliar integers

! Allocation arrays to stock matrix morse
! Array number of last term of every line (mors1)
   subroutine initialize_poisson_solver(self, mesh, ntypfr, potfr)

      type(sll_t_poisson_2d_triangular), intent(out)         :: self
      type(sll_t_triangular_mesh_2d), intent(in), target :: mesh
      sll_int32, intent(in)         :: ntypfr(:)
      sll_real64, intent(in)         :: potfr(:)

      sll_real64, allocatable     :: tmp1(:)
      character(len=*), parameter :: self_sub_name = 'read_data_solver'
      character(len=128)          :: err_msg

      sll_int32 :: nref, nn, ndir
      sll_int32 :: i
      sll_int32 :: ierr

      if (mesh%analyzed) then
         self%mesh => mesh
      else
         err_msg = "Call sll_s_analyze_triangular_mesh before initialize_poisson_solver."
         SLL_ERROR(self_sub_name, err_msg)
      end if

      self%ntypfr = ntypfr
      self%potfr = potfr

      allocate (self%sv1(mesh%nbtcot)); self%sv1 = 0.0_f64
      allocate (self%sv2(mesh%nbtcot)); self%sv2 = 0.0_f64
      allocate (self%vtantx(mesh%nbtcot)); self%vtantx = 0.0_f64
      allocate (self%vtanty(mesh%nbtcot)); self%vtanty = 0.0_f64

      allocate (self%mors1(mesh%num_nodes + 1)); self%mors1 = 0

! Array witth number of every line term (mors2)

      allocate (self%mors2(12*mesh%num_nodes)); self%mors2 = 0

! mors1,mors2.

      call morse(mesh%npoel1, &
                 mesh%npoel2, &
                 mesh%nodes, &
                 mesh%num_triangles, &
                 mesh%num_nodes, &
                 self%mors1, &
                 self%mors2)

      allocate (self%iprof(mesh%num_nodes + 1)); self%iprof = 0

      call profil(mesh%nodes, &
                  mesh%num_nodes, &
                  mesh%npoel1, &
                  mesh%npoel2, &
                  self%iprof)

!======================================================================
!--- 2.0 --- POISSON finite element -----------------
!======================================================================

!matrix de masse diagonalized
      allocate (self%amass(mesh%num_nodes)); self%amass = 0.0_f64

!matrix "grad-grad" stocked profil form.
      allocate (self%grgr(self%iprof(mesh%num_nodes + 1))); self%grgr = 0.0_f64

!gradx grady
      allocate (self%gradx(self%mors1(mesh%num_nodes + 1))); self%gradx = 0.0_f64
      allocate (self%grady(self%mors1(mesh%num_nodes + 1))); self%grady = 0.0_f64

!--- Array boundaries Dirichlet -------------------------

      ndir = 0
      do nn = 1, mesh%num_nodes
         nref = mesh%refs(nn)
         if (nref > 0) then
            if (self%ntypfr(nref) == 1) then
               ndir = ndir + 1
            end if
         end if
      end do

      allocate (self%ifron(ndir)); self%ifron = 0

      ndir = 0
      do nn = 1, mesh%num_nodes
         nref = mesh%refs(nn)
         if (nref > 0) then
            if (self%ntypfr(nref) == 1) then
               ndir = ndir + 1
               self%ifron(ndir) = nn
            end if
         end if
      end do

      self%ndiric = ndir

!--- Calcul des matrix ----------------------------------------------

      call poismc(self)

!Calcul de la matrix B tel que B*Bt = A dans le cas Cholesky

      allocate (tmp1(self%iprof(mesh%num_nodes + 1))); tmp1 = 0.0_f64

!write(*,"(//5x,a)")" *** Appel Choleski pour Poisson ***  "
      call sll_s_choles(self%iprof, self%grgr, tmp1)

      do i = 1, self%iprof(mesh%num_nodes + 1)
         self%grgr(i) = tmp1(i)
      end do

      deallocate (tmp1)

      SLL_CLEAR_ALLOCATE(self%vnx(1:self%mesh%num_nodes), ierr)
      SLL_CLEAR_ALLOCATE(self%vny(1:self%mesh%num_nodes), ierr)
      SLL_ALLOCATE(self%naux(1:self%mesh%num_nodes), ierr)
      self%naux = .false.

   end subroutine initialize_poisson_solver

!======================================================================

! Function: morse
! Calculation of the arrays containing the address of the last term of
! each line of the "walrus" matrix and the arrays containing the
! numbers of the terms of these matrixes. The diagonal term is
! last place in each line.
!
! npoel1 - location in npoel2 of the last element relative to each
! node with npoel1(1)=0 and npoel1(i+1) relative to node i
! npoel2 - arrays numbers of elements having a common vertex
! ntri - triangles vertex numbers
! nbt - number of mesh triangles
! nbs - number of nodes
!
! mors1 - arrays addresses of the last terms of each line
! the matrix with the convention:
! mors1(1)=0 and mors1(i+1) aii address
! mors2 - arrays of the numbers of the "morse" matrix terms

   subroutine morse(npoel1, npoel2, ntri, nbt, nbs, mors1, mors2)

      sll_int32, intent(in)  :: nbt
      sll_int32, intent(in)  :: nbs
      sll_int32, dimension(:), intent(in)  :: npoel1
      sll_int32, dimension(:), intent(in)  :: npoel2
      sll_int32, dimension(3, nbt), intent(in)  :: ntri
      sll_int32, dimension(:), intent(out) :: mors1
      sll_int32, dimension(:), intent(out) :: mors2
      sll_int32, dimension(20)                 :: ilign

      sll_int32 :: l, itest1, itest2, js1, js2, is1, is2, is3, numel
      sll_int32 :: iel, nlign, nel, is, im, k

      im = 0
      k = 0

      mors1(1) = 0

      do is = 1, nbs  !Boucle sur les nodes

         !elements 'is' node

         nel = npoel1(is + 1) - npoel1(is)
         nlign = 0

         !Loop over this elements

         do iel = 1, nel

            k = k + 1

            numel = npoel2(k)
            is1 = ntri(1, numel); is2 = ntri(2, numel); is3 = ntri(3, numel)
            if (is1 == is) then
               js1 = is2; js2 = is3
            end if
            if (is2 == is) then
               js1 = is3; js2 = is1
            end if
            if (is3 == is) then
               js1 = is1; js2 = is2
            end if

            !We see if the 2 nodes other than is of the current element
            !have already been taken into account in all the interacting nodes
            !with is.

            itest1 = 0; itest2 = 0
            if (nlign .NE. 0) then
               do l = 1, nlign
                  if (js1 == ilign(l)) then
                     itest1 = 1
                  end if
                  if (js2 == ilign(l)) then
                     itest2 = 1
                  end if
               end do
            end if

            if (itest1 == 0) then
               nlign = nlign + 1
               ilign(nlign) = js1
            end if
            if (itest2 == 0) then
               nlign = nlign + 1
               ilign(nlign) = js2
            end if

         end do

         !Definition of the address of the last term of the line
         mors1(is + 1) = mors1(is) + nlign + 1

         !Filling the jaw arrays2 with the term numbers
         !of the is line.

         if (nlign .NE. 0) then
            do l = 1, nlign
               im = im + 1
               mors2(im) = ilign(l)
            end do
         end if
         im = im + 1
         mors2(im) = is
      end do

   end subroutine morse

! ======================================================================

!Function: poismc
!
!   Compute all matrices for cartesian Poisson
!
! Input:
!
! m      - super-arrays
! coor   - coord nodes
! refs   - indice if node in ref boundary
! ifron  - arrays of nodes Dirichlet
! mors1  - arrays of number terms per line matrix morse symetric
! mors2  - arrays of number terms every line matrix morse symetric
! ntri   - numbers triangles nodes
! area   - area triangles
! iprof  - profile matrix grad-grad
!
! noefnt - nodes interns Dirichlet
!
! Output:
!
! agrgr  - matrix operator "grad-grad" "morse"
! aliss  - matrix fitting
! amass  - matrix mass diag
! amclt  - matrix mass diag absorbing bc
! aroro  - matrix operator "rot-rot" form "morse"
! d1dx,  - deriv  functions basis every triangle
!..d3dy
! gradx  - matrix operator gradx
! grady  - matrix operator grady
! rotx   - matrix operator rotx
! roty   - matrix operator roty
! area   - every element
!
!Auteur: Puertolas - Version 1.0  Octobre  1992
   subroutine poismc(self)

      type(sll_t_poisson_2d_triangular), intent(inout) :: self
      sll_real64 :: amloc(3), aggloc(9), grxloc(9), gryloc(9)
      sll_real64 :: dntx1, dntx2, dntx3, dnty1, dnty2, dnty3
      sll_real64 :: x1t, x2t, x3t, y1t, y2t, y3t, coef
      sll_int32 :: is1t, is2t, is3t, iel
      sll_int32 :: is, j

!Boucle sur les elements.

      do iel = 1, self%mesh%num_triangles

         !coefficients geometry

         is1t = self%mesh%nodes(1, iel)
         is2t = self%mesh%nodes(2, iel)
         is3t = self%mesh%nodes(3, iel)

         x1t = self%mesh%coord(1, is1t)
         x2t = self%mesh%coord(1, is2t)
         x3t = self%mesh%coord(1, is3t)

         y1t = self%mesh%coord(2, is1t)
         y2t = self%mesh%coord(2, is2t)
         y3t = self%mesh%coord(2, is3t)

         dntx1 = y2t - y3t
         dntx2 = y3t - y1t
         dntx3 = y1t - y2t

         dnty1 = x3t - x2t
         dnty2 = x1t - x3t
         dnty3 = x2t - x1t

         !Contribution matrix mass

         amloc(1) = self%mesh%area(iel)/3.
         amloc(2) = self%mesh%area(iel)/3.
         amloc(3) = self%mesh%area(iel)/3.

         !Assembling

         call asbld(amloc, is1t, is2t, is3t, self%amass)

         !Contribution matrix grad-grad

         coef = 1./(4.*self%mesh%area(iel))

         aggloc(1) = (dntx1**2 + dnty1**2)*coef
         aggloc(2) = (dntx1*dntx2 + dnty1*dnty2)*coef
         aggloc(3) = (dntx1*dntx3 + dnty1*dnty3)*coef
         aggloc(4) = (dntx2**2 + dnty2**2)*coef
         aggloc(5) = (dntx2*dntx3 + dnty2*dnty3)*coef
         aggloc(6) = (dntx3**2 + dnty3**2)*coef

         !Contribution matrix gradx grady:

         grxloc(1) = -dntx1/6.; gryloc(1) = -dnty1/6.
         grxloc(2) = -dntx2/6.; gryloc(2) = -dnty2/6.
         grxloc(3) = -dntx3/6.; gryloc(3) = -dnty3/6.
         grxloc(4) = -dntx1/6.; gryloc(4) = -dnty1/6.
         grxloc(5) = -dntx2/6.; gryloc(5) = -dnty2/6.
         grxloc(6) = -dntx3/6.; gryloc(6) = -dnty3/6.
         grxloc(7) = -dntx1/6.; gryloc(7) = -dnty1/6.
         grxloc(8) = -dntx2/6.; gryloc(8) = -dnty2/6.
         grxloc(9) = -dntx3/6.; gryloc(9) = -dnty3/6.

         !Assembling

         call asblp(self%iprof, aggloc, is1t, is2t, is3t, self%grgr)

         call asblm2(grxloc, &
                     gryloc, &
                     self%mors1, &
                     self%mors2, &
                     is1t, &
                     is2t, &
                     is3t, &
                     self%gradx, &
                     self%grady)

      end do

! ======================================================================
! ... Dirichlet

      do j = 1, self%ndiric
         is = self%ifron(j)
         self%grgr(self%iprof(is + 1)) = grandx
      end do

   end subroutine poismc

!Function: poissn
!potential and electric field
!Poisson en cartesien
!
! ex     - e1  nodes
! ey     - e2  nodes
! rho    - density charge nodes
! phi    - potential nodes
! ifron  - nodes boundary Dirichlet
! noefnt - nodes internal boundaries Dirichlet
! irffnt - reference of this nodes
! grgr   - matrix operator "grad-grad"
! amass  - matrix mass diag
! gra1   - matrix operator "grad1"
! gra2   - matrix operator "grad2"
! mors1  - arrays matrix "morse"
! mors2  - arrays matrix "morse"
! iprof  - matrix profil
!
!
! grgrdd - product matrix grgr vector p (gradient conjugue)
! dird   - direction method gradient conjugue
! res    - residual method gradient conjugue
! sdmb   - second member system potential
! sdmb12 - second member system electric field
! precd  - arrays local utilise gradient conjugue preconditionned
! nbs    - number nodes mesh
! nmxfr  - number boundaries ref
! errpoi - precision
! nitgc  - number max iterations gradient conjugue
!
   subroutine poissn(self, rho, phi, ex, ey)

      type(sll_t_poisson_2d_triangular)                :: self
      sll_real64, dimension(:), intent(in)            :: rho
      sll_real64, dimension(:), intent(out)           :: phi
      sll_real64, dimension(:), intent(out), optional :: ex
      sll_real64, dimension(:), intent(out), optional :: ey

      sll_real64, dimension(:), allocatable           :: sdmb
      sll_real64, dimension(:), allocatable           :: sdmb12

      sll_int32 :: i, is, nref, ierr

! ----------- CALCUL DU POTENTIEL  -------------------------------------
!
!... RHS ...

      SLL_ALLOCATE(sdmb(self%mesh%num_nodes), ierr)

!!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(is)
!!$OMP DO SCHEDULE(RUNTIME)
      do is = 1, self%mesh%num_nodes
         sdmb(is) = self%amass(is)*rho(is)/self%eps0
      end do
!!$OMP END DO
!!$OMP END PARALLEL

!... Dirichlet

      do is = 1, self%ndiric
         nref = self%mesh%refs(self%ifron(is))
         sdmb(self%ifron(is)) = self%potfr(nref)*grandx
      end do

      call sll_s_desrem(self%iprof, self%grgr, sdmb, self%mesh%num_nodes, phi)

!*** Ex  Ey:

      if (present(ex) .and. present(ey)) then

         SLL_ALLOCATE(sdmb12(self%mesh%num_nodes), ierr)
         !*** RHS for Ex:

         call m1p(self%gradx, self%mors1, self%mors2, phi, self%mesh%num_nodes, sdmb12)

         do i = 1, self%mesh%num_nodes
            ex(i) = sdmb12(i)/self%amass(i)
         end do

         !*** RHS for Ey:

         call m1p(self%grady, self%mors1, self%mors2, phi, self%mesh%num_nodes, sdmb12)

         do i = 1, self%mesh%num_nodes
            ey(i) = sdmb12(i)/self%amass(i)
         end do

      end if

   end subroutine poissn

!Function: poifrc
!Correction on boundaries
!
!  - E.tau = 0  boundaries Dirichlets
!  - E.nu =  0  boundaries Neumann
!
!Input:
!
!  ksofro  - numbers edge nodes
!  krefro  - numbers edge reference
!  vnofro  - vector normal (inside)
!  vnx,vny - arrays normal vectors nodes
!
!Array local:
!
! naux   - arrays for nodes on boundary
!
   subroutine poifrc(self, ex, ey)

      type(sll_t_poisson_2d_triangular)   :: self
      sll_real64 :: ex(:), ey(:)
      sll_real64 :: pscal, xnor
      sll_int32 :: is1, is2, ict
      sll_int32 :: i

!!$ ======================================================================
!!$ ... E.tau = 0 boundaries Dirichlet

!!$ ... loop over edges boundaries for normals
!!$ ... nodes "Dirichlet"
      self%vnx = 0.0_f64; self%vny = 0.0_f64; self%naux = .false.

      do ict = 1, self%mesh%nctfrt

         if (self%ntypfr(self%mesh%krefro(ict)) == 1) then

            is1 = self%mesh%ksofro(1, ict)
            is2 = self%mesh%ksofro(2, ict)

            self%vnx(is1) = self%vnx(is1) + self%mesh%vnofro(1, ict)
            self%vny(is1) = self%vny(is1) + self%mesh%vnofro(2, ict)
            self%vnx(is2) = self%vnx(is2) + self%mesh%vnofro(1, ict)
            self%vny(is2) = self%vny(is2) + self%mesh%vnofro(2, ict)

            self%naux(is1) = .true.
            self%naux(is2) = .true.

         end if

      end do

!... E.tau=0

      do i = 1, self%mesh%num_nodes

         if (self%naux(i)) then

            xnor = SQRT(self%vnx(i)**2 + self%vny(i)**2)
            if (xnor > self%mesh%petitl) then
               self%vnx(i) = self%vnx(i)/xnor
               self%vny(i) = self%vny(i)/xnor

               pscal = self%vnx(i)*ex(i) + self%vny(i)*ey(i)
               ex(i) = self%vnx(i)*pscal
               ey(i) = self%vny(i)*pscal
            end if

         end if

      end do

!======================================================================
! ... E.nu = 0 boundaries Neumann
!
      self%vnx = 0.0_f64; self%vny = 0.0_f64; self%naux = .false.

! ... Loop over edges boundaries nodes "Neumann"

      do ict = 1, self%mesh%nctfrt

         if (self%ntypfr(self%mesh%krefro(ict)) == 3) then

            is1 = self%mesh%ksofro(1, ict)
            is2 = self%mesh%ksofro(2, ict)

            self%vnx(is1) = self%vnx(is1) + self%mesh%vnofro(1, ict)
            self%vny(is1) = self%vny(is1) + self%mesh%vnofro(2, ict)
            self%vnx(is2) = self%vnx(is2) + self%mesh%vnofro(1, ict)
            self%vny(is2) = self%vny(is2) + self%mesh%vnofro(2, ict)

            self%naux(is1) = .true.
            self%naux(is2) = .true.

         end if

      end do

!... E.nu=0

      do i = 1, self%mesh%num_nodes

         if (self%naux(i)) then

            xnor = SQRT(self%vnx(i)**2 + self%vny(i)**2)
            if (xnor > self%mesh%petitl) then
               self%vnx(i) = self%vnx(i)/xnor
               self%vny(i) = self%vny(i)/xnor
               pscal = self%vnx(i)*ex(i) + self%vny(i)*ey(i)
               ex(i) = ex(i) - self%vnx(i)*pscal
               ey(i) = ey(i) - self%vny(i)*pscal
            end if

         end if

      end do

   end subroutine poifrc

!Function: profil
!numbers elements on diagonal matrix profil associate
!
! Input:
! npoel1 - index in npoel2 of last element of every node
!          with npoel1(1)=0 and npoel1(i+1) link to node i
! npoel2 - arrays numbers elements with common node
! ntri   - numbers nodes triangles
! nbt    - number triangles
! noe    - number nodes
!
! Ouput:
! iprof  - arrays numbers terms diagonal with
!        - iprof(i+1) term diag line i
!        - iprof(1)=0
!
!Auteur:
! J. Segre - Juillet 89
   subroutine profil(ntri, nbs, npoel1, npoel2, iprof)

      sll_int32, dimension(:) :: iprof
      sll_int32, dimension(:), intent(in) :: npoel1, npoel2
      sll_int32, dimension(:, :), intent(in) :: ntri
      sll_int32 :: in, k, is1, is2, is3, numel, ind, iel, nel
      sll_int32, intent(in) :: nbs

!************************************************************************
!*
!* Compute length of line i matrix profil :
!* for every node i we look for nodes j linked;
!* if node j not found, length line j is set to
!* j-i+1
!
!* loop over nodes

      k = 0
      do in = 1, nbs

         !* number elements with in node
         nel = npoel1(in + 1) - npoel1(in)

         !* loop over these elements

         do iel = 1, nel

            k = k + 1
            !* number element
            numel = npoel2(k)
            is1 = ntri(1, numel)
            ind = is1 + 1
            if (iprof(ind) .eq. 0) then
               iprof(ind) = is1 - in + 1
            end if
            is2 = ntri(2, numel)
            ind = is2 + 1
            if (iprof(ind) .eq. 0) then
               iprof(ind) = is2 - in + 1
            end if
            is3 = ntri(3, numel)
            ind = is3 + 1
            if (iprof(ind) .eq. 0) then
               iprof(ind) = is3 - in + 1
            end if

         end do

      end do

!* compute position terms diags matrix
!* profil (sum numbers elements of lines
!* before and current line ).

      do ind = 3, nbs + 1
         iprof(ind) = iprof(ind - 1) + iprof(ind)
      end do

   end subroutine profil

!Function: asbld
!  Assemble matrix element in matrix global if diagonal
!
!Parametres d'entree:
!
!  aele         -    Matrix diag
!                    (3 terms per element triangular)
!  i1,i2,i3     -    numbers nodes element
!
!Parametre resultat:
!
!  xmass    -   matrix global diagonal
!
!Auteur:
!      J. Segre - Version 1.0  Juillet 1989
   subroutine asbld(aele, i1, i2, i3, xmass)

      sll_real64, dimension(:) :: aele, xmass
      sll_int32 :: i1, i2, i3

      xmass(i1) = xmass(i1) + aele(1)
      xmass(i2) = xmass(i2) + aele(2)
      xmass(i3) = xmass(i3) + aele(3)

   end subroutine

!Function: asblm2
!
!          assemble 3 matrix element
!          3 matrix global stocked
!          form "morse" non symetric
!
!      Input:
!          aele1,aele2    -  Matrices elements
!          mors1          -  Array number of terms per
!                            line matrix morse
!          mors2          -  Array numbers terms
!                            of every line matrix morse
!          i1,i2,i3       -  Num nodes element
!
!      Output:
!          a1,a2          -  Matrix global stocked
!                            format "morse"
!
!Auteur:
!J. Segre - Version 1.0  Aout 1989
   subroutine asblm2(aele1, aele2, mors1, mors2, i1, i2, i3, a1, a2)

      sll_int32, intent(in) :: i1, i2, i3
      sll_int32, dimension(:), intent(in) :: mors1, mors2
      sll_real64, dimension(:), intent(in)  :: aele1, aele2
      sll_real64, dimension(:), intent(out) :: a1, a2
      sll_int32 :: j, j1, j2, ind1, ind2, ind3

! --- 1.1 --- terms diag ---------------------------

      ind1 = mors1(i1 + 1)
      a1(ind1) = a1(ind1) + aele1(1)
      a2(ind1) = a2(ind1) + aele2(1)

      ind2 = mors1(i2 + 1)
      a1(ind2) = a1(ind2) + aele1(5)
      a2(ind2) = a2(ind2) + aele2(5)

      ind3 = mors1(i3 + 1)
      a1(ind3) = a1(ind3) + aele1(9)
      a2(ind3) = a2(ind3) + aele2(9)

! --- 1.2 --- other terms ------------------------------

      j2 = ind1 - 1
      j1 = mors1(i1) + 1

      if (j2 >= j1) then
         do j = j1, j2
            if (i2 == mors2(j)) then
               a1(j) = a1(j) + aele1(2)
               a2(j) = a2(j) + aele2(2)
            end if
            if (i3 == mors2(j)) then
               a1(j) = a1(j) + aele1(3)
               a2(j) = a2(j) + aele2(3)
            end if
         end do
      end if

      j2 = ind2 - 1
      j1 = mors1(i2) + 1

      if (j2 >= j1) then
         do j = j1, j2
            if (i1 == mors2(j)) then
               a1(j) = a1(j) + aele1(4)
               a2(j) = a2(j) + aele2(4)
            end if
            if (i3 == mors2(j)) then
               a1(j) = a1(j) + aele1(6)
               a2(j) = a2(j) + aele2(6)
            end if
         end do
      end if

      j2 = ind3 - 1
      j1 = mors1(i3) + 1

      if (j2 >= j1) then
         do j = j1, j2
            if (i1 == mors2(j)) then
               a1(j) = a1(j) + aele1(7)
               a2(j) = a2(j) + aele2(7)
            end if
            if (i2 == mors2(j)) then
               a1(j) = a1(j) + aele1(8)
               a2(j) = a2(j) + aele2(8)
            end if
         end do
      end if

   end subroutine

!----------------------------------------------------------------------

!Function: asblp
!
!      assemble matrix element
!      in matrix global symetric
!      format "profil"
!
!Input:
!
!          aele         -  Matrix element (6 - 10 terms
!                          element triangle or quad)
!          iprof        -  Description matrix
!          i1,i2,i3     -  Num nodes element
!
!Output:
!
!          xmass        -  matrix global stockee format "profil"
!
!Auteur:
!       J. Segre - Version 1.0  Aout 1989
   subroutine asblp(iprof, aele, i1, i2, i3, xmass)

      sll_int32, dimension(:) :: iprof
      sll_real64, dimension(:) :: aele(*), xmass(*)
      sll_int32 :: i1, i2, i3, idiag1, idiag2, idiag3, ind

!--- 1.1 --- terms diag ---------------------------

      idiag1 = iprof(i1 + 1)
      idiag2 = iprof(i2 + 1)
      idiag3 = iprof(i3 + 1)

      xmass(idiag1) = xmass(idiag1) + aele(1)
      xmass(idiag2) = xmass(idiag2) + aele(4)
      xmass(idiag3) = xmass(idiag3) + aele(6)

!--- 1.2 --- others terms ------------------------------
!           (save only aij if i>j)

      if (i1 < i2) then
         ind = idiag2 + i1 - i2
      else
         ind = idiag1 + i2 - i1
      end if

      xmass(ind) = xmass(ind) + aele(2)

      if (i1 < i3) then
         ind = idiag3 + i1 - i3
      else
         ind = idiag1 + i3 - i1
      end if

      xmass(ind) = xmass(ind) + aele(3)

      if (i2 < i3) then
         ind = idiag3 + i2 - i3
      else
         ind = idiag2 + i3 - i2
      end if

      xmass(ind) = xmass(ind) + aele(5)

!----------------------------------------------------------------------

   end subroutine asblp

!Function: m1p
!
!     yvect =  xmors.xvect
!     xvect,yvect  vectors
!
!Input:
!          xmors        -     matrix  "morse"
!          xvect        -     vector
!          mors1,mors2  -     arrays matrix "morse"
!          nlign        -     number lines matrices
!
!Output:
!          yvect    -   vector result
!
!Auteur:
!     J. Segre - Version 1.0  Decembre 1989
   subroutine m1p(xmors, mors1, mors2, xvect, nlign, yvect)

      sll_real64    :: xmors(*), xvect(*), yvect(*)
      sll_int32 :: mors1(*), mors2(*)
      sll_int32 :: il, nlign, noeui

      do il = 1, nlign

         noeui = mors1(il + 1) - mors1(il)

         select case (noeui)
         case (6)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5))

         case (5)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4))

         case (7)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5)) + &
                        xmors(mors1(il + 1) - 6)*xvect(mors2(mors1(il + 1) - 6))

         case (4)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3))

         case (8)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5)) + &
                        xmors(mors1(il + 1) - 6)*xvect(mors2(mors1(il + 1) - 6)) + &
                        xmors(mors1(il + 1) - 7)*xvect(mors2(mors1(il + 1) - 7))

         case (3)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2))

         case (9)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5)) + &
                        xmors(mors1(il + 1) - 6)*xvect(mors2(mors1(il + 1) - 6)) + &
                        xmors(mors1(il + 1) - 7)*xvect(mors2(mors1(il + 1) - 7)) + &
                        xmors(mors1(il + 1) - 8)*xvect(mors2(mors1(il + 1) - 8))

         case (10)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5)) + &
                        xmors(mors1(il + 1) - 6)*xvect(mors2(mors1(il + 1) - 6)) + &
                        xmors(mors1(il + 1) - 7)*xvect(mors2(mors1(il + 1) - 7)) + &
                        xmors(mors1(il + 1) - 8)*xvect(mors2(mors1(il + 1) - 8)) + &
                        xmors(mors1(il + 1) - 9)*xvect(mors2(mors1(il + 1) - 9))

         case (11)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5)) + &
                        xmors(mors1(il + 1) - 6)*xvect(mors2(mors1(il + 1) - 6)) + &
                        xmors(mors1(il + 1) - 7)*xvect(mors2(mors1(il + 1) - 7)) + &
                        xmors(mors1(il + 1) - 8)*xvect(mors2(mors1(il + 1) - 8)) + &
                        xmors(mors1(il + 1) - 9)*xvect(mors2(mors1(il + 1) - 9)) + &
                        xmors(mors1(il + 1) - 10)*xvect(mors2(mors1(il + 1) - 10))

         case (12)
            yvect(il) = xmors(mors1(il + 1))*xvect(mors2(mors1(il + 1))) + &
                        xmors(mors1(il + 1) - 1)*xvect(mors2(mors1(il + 1) - 1)) + &
                        xmors(mors1(il + 1) - 2)*xvect(mors2(mors1(il + 1) - 2)) + &
                        xmors(mors1(il + 1) - 3)*xvect(mors2(mors1(il + 1) - 3)) + &
                        xmors(mors1(il + 1) - 4)*xvect(mors2(mors1(il + 1) - 4)) + &
                        xmors(mors1(il + 1) - 5)*xvect(mors2(mors1(il + 1) - 5)) + &
                        xmors(mors1(il + 1) - 6)*xvect(mors2(mors1(il + 1) - 6)) + &
                        xmors(mors1(il + 1) - 7)*xvect(mors2(mors1(il + 1) - 7)) + &
                        xmors(mors1(il + 1) - 8)*xvect(mors2(mors1(il + 1) - 8)) + &
                        xmors(mors1(il + 1) - 9)*xvect(mors2(mors1(il + 1) - 9)) + &
                        xmors(mors1(il + 1) - 10)*xvect(mors2(mors1(il + 1) - 10)) + &
                        xmors(mors1(il + 1) - 11)*xvect(mors2(mors1(il + 1) - 11))

         end select

      end do

   end subroutine m1p

!  Compute Ex et Ey from potential
   subroutine poliss(self, phi, ex, ey)

      type(sll_t_poisson_2d_triangular), intent(inout) :: self
      sll_real64, dimension(:), intent(in)    :: phi
      sll_real64, dimension(:), intent(out)   :: ex
      sll_real64, dimension(:), intent(out)   :: ey

      sll_int32 :: is, ic, nbc, iac
      LOGICAL :: lerr

      lerr = .FALSE.

      do ic = 1, self%mesh%nbtcot
         self%vtanty(ic) = (phi(self%mesh%nuvac(1, ic)) - phi(self%mesh%nuvac(2, ic)))/self%mesh%xlcod(ic)
      end do

      do ic = 1, self%mesh%nbtcot
         self%vtantx(ic) = self%vtanty(ic)*self%mesh%vtaux(ic)
      end do

      do ic = 1, self%mesh%nbtcot
         self%vtanty(ic) = self%vtanty(ic)*self%mesh%vtauy(ic)
      end do

      do is = 1, self%mesh%num_nodes

         iac = self%mesh%nbcov(is) + 1
         nbc = self%mesh%nbcov(is + 1) - self%mesh%nbcov(is)

         if (nbc == 6) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5))

         else if (nbc == 5) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4))

         else if (nbc == 7) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5)) &
                           + self%vtantx(self%mesh%nugcv(iac + 6))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5)) &
                           + self%vtanty(self%mesh%nugcv(iac + 6))

         else if (nbc == 4) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3))

         else if (nbc == 8) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5)) &
                           + self%vtantx(self%mesh%nugcv(iac + 6)) &
                           + self%vtantx(self%mesh%nugcv(iac + 7))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5)) &
                           + self%vtanty(self%mesh%nugcv(iac + 6)) &
                           + self%vtanty(self%mesh%nugcv(iac + 7))

         else if (nbc == 3) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2))

         else if (nbc == 9) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5)) &
                           + self%vtantx(self%mesh%nugcv(iac + 6)) &
                           + self%vtantx(self%mesh%nugcv(iac + 7)) &
                           + self%vtantx(self%mesh%nugcv(iac + 8))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5)) &
                           + self%vtanty(self%mesh%nugcv(iac + 6)) &
                           + self%vtanty(self%mesh%nugcv(iac + 7)) &
                           + self%vtanty(self%mesh%nugcv(iac + 8))

         else if (nbc == 2) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1))

         else if (nbc == 10) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5)) &
                           + self%vtantx(self%mesh%nugcv(iac + 6)) &
                           + self%vtantx(self%mesh%nugcv(iac + 7)) &
                           + self%vtantx(self%mesh%nugcv(iac + 8)) &
                           + self%vtantx(self%mesh%nugcv(iac + 9))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5)) &
                           + self%vtanty(self%mesh%nugcv(iac + 6)) &
                           + self%vtanty(self%mesh%nugcv(iac + 7)) &
                           + self%vtanty(self%mesh%nugcv(iac + 8)) &
                           + self%vtanty(self%mesh%nugcv(iac + 9))

         else if (nbc == 11) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5)) &
                           + self%vtantx(self%mesh%nugcv(iac + 6)) &
                           + self%vtantx(self%mesh%nugcv(iac + 7)) &
                           + self%vtantx(self%mesh%nugcv(iac + 8)) &
                           + self%vtantx(self%mesh%nugcv(iac + 9)) &
                           + self%vtantx(self%mesh%nugcv(iac + 10))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5)) &
                           + self%vtanty(self%mesh%nugcv(iac + 6)) &
                           + self%vtanty(self%mesh%nugcv(iac + 7)) &
                           + self%vtanty(self%mesh%nugcv(iac + 8)) &
                           + self%vtanty(self%mesh%nugcv(iac + 9)) &
                           + self%vtanty(self%mesh%nugcv(iac + 10))

         else if (nbc == 12) then

            self%sv1(is) = self%vtantx(self%mesh%nugcv(iac)) &
                           + self%vtantx(self%mesh%nugcv(iac + 1)) &
                           + self%vtantx(self%mesh%nugcv(iac + 2)) &
                           + self%vtantx(self%mesh%nugcv(iac + 3)) &
                           + self%vtantx(self%mesh%nugcv(iac + 4)) &
                           + self%vtantx(self%mesh%nugcv(iac + 5)) &
                           + self%vtantx(self%mesh%nugcv(iac + 6)) &
                           + self%vtantx(self%mesh%nugcv(iac + 7)) &
                           + self%vtantx(self%mesh%nugcv(iac + 8)) &
                           + self%vtantx(self%mesh%nugcv(iac + 9)) &
                           + self%vtantx(self%mesh%nugcv(iac + 10)) &
                           + self%vtantx(self%mesh%nugcv(iac + 11))

            self%sv2(is) = self%vtanty(self%mesh%nugcv(iac)) &
                           + self%vtanty(self%mesh%nugcv(iac + 1)) &
                           + self%vtanty(self%mesh%nugcv(iac + 2)) &
                           + self%vtanty(self%mesh%nugcv(iac + 3)) &
                           + self%vtanty(self%mesh%nugcv(iac + 4)) &
                           + self%vtanty(self%mesh%nugcv(iac + 5)) &
                           + self%vtanty(self%mesh%nugcv(iac + 6)) &
                           + self%vtanty(self%mesh%nugcv(iac + 7)) &
                           + self%vtanty(self%mesh%nugcv(iac + 8)) &
                           + self%vtanty(self%mesh%nugcv(iac + 9)) &
                           + self%vtanty(self%mesh%nugcv(iac + 10)) &
                           + self%vtanty(self%mesh%nugcv(iac + 11))
         else

            lerr = .TRUE.

         end if

      end do

      if (lerr) then
         write (6, 900)
         stop
      end if

! --- 3.0 --- Solve systems 2*2 --------------------

      do is = 1, self%mesh%num_nodes
         ex(is) = self%mesh%xmal2(is)*self%sv1(is) - self%mesh%xmal3(is)*self%sv2(is)
         ey(is) = self%mesh%xmal1(is)*self%sv2(is) - self%mesh%xmal3(is)*self%sv1(is)
      end do

! --- 9.0 --- Formats --------------------------------------------------

900   format(//10x, 'Error dans POLISS' &
              /10x, 'On a trouve more de 12 edges')

   end subroutine poliss

end module sll_m_poisson_2d_tri

