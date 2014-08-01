program test_hex_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use hex_mesh
  use sll_hermite_interpolation_2d_module

  implicit none

  type(deriv), pointer        :: deriv 
  type(hex_mesh_2d), pointer  :: mesh


  sll_int32    :: num_cells
  sll_int32    :: i
  sll_int32    :: nloops
  sll_int32    :: ierr
  ! initial distribution
  sll_real64   :: gauss_x2
  sll_real64   :: gauss_x1
  sll_real64   :: gauss_sig
  sll_real64,dimension(:),allocatable :: x1
  sll_real64,dimension(:),allocatable :: x2
  sll_real64,dimension(:),allocatable :: f_init
  sll_real64,dimension(:),allocatable :: chi1
  sll_real64,dimension(:),allocatable :: chi2
  ! distribution at time n
  sll_real64,dimension(:),allocatable :: x1_char
  sll_real64,dimension(:),allocatable :: x2_char
  sll_real64,dimension(:),allocatable :: f_tn
  ! distribution at end time
  sll_real64,dimension(:),allocatable :: f_fin
  sll_int32    :: where_error
  sll_real64   :: diff_error
  sll_real64   :: norm2_error
  ! advection
  sll_int32    :: which_advec
  sll_real64   :: advec
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  ! timer variables
  sll_real64   :: tcpu, t_init, t_end
  !others
  sll_real64   :: x2_basis
  sll_real64   :: x1_basis
  sll_real64   :: x1_temp
  character(len = 50) :: filename
  character(len = 50) :: filename2
  character(len = 4)  :: filenum
  sll_real64   :: p = 6 !-> degree of the approximation for the derivative 



  do num_cells = 75,75,1 ! -> boucle sur la taille du maillage  

     !   t_init = getRealTimer()
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Mesh initialization   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     mesh => new_hex_mesh_2d(num_cells, 0._f64, 0._f64, &
                                !0.5_f64, -sqrt(3.)/2._f64, &
                                !0.5_f64,  sqrt(3.)/2._f64, &
                                !1.0_f64,  0._f64, &
          radius = 1._f64)

     call sll_display(mesh) !  -> décrit sur le terminal le mesh

     print*,"num_pts : ", mesh%num_pts_tot

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !             allocation
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     SLL_ALLOCATE(f_init(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(f_tn(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(f_fin(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(x1(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(x2(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(x1_char(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(x2_char(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(chi1(mesh%num_pts_tot),ierr)
     SLL_ALLOCATE(chi2(mesh%num_pts_tot),ierr)

     SLL_ALLOCATE(deriv(mesh%num_pts_tot),ierr)

     
     a    = radius / real( num_cells )
     aire = a**2/sqrt(3_f64)

     ! Distribution initialization

     gauss_x1  = -0.25_f64
     gauss_x2  = -0.25_f64
     gauss_sig = 0.05_f64

     do i=0, mesh%num_pts_tot-1
        x1(i+1) = from_global_x1(mesh, i)
        x2(i+1) = from_global_x2(mesh, i)
        f_init(i+1) = 1._f64*exp(-0.5_f64*((x1(i+1)-gauss_x1)**2 / gauss_sig**2 &
             + (x2(i+1)-gauss_x2)**2 / gauss_sig**2))
        if (exponent(f_init(i+1)) .lt. -17) then
           f_init(i+1) = 0._f64
        end if
        f_tn(i+1) = f_init(i+1)
     end do

     call write_field_hex_mesh(mesh, f_init, "init_dist.txt")

     ! Advection initialization
     which_advec = 1
     advec = 0.025_f64!5_f64
     tmax  = 2.5_f64
     dt    = 0.025_f64
     t     = 0._f64

     ! Computing characteristics
     if (which_advec .eq. 0) then
        ! linear advection
        x1_char(:) = x1(:) - advec*dt
        x2_char(:) = x2(:) - advec*dt
     else
        ! Circular advection
        x1_char(1) = 0._f64
        x2_char(1) = 0._f64
        x1_char(2:) = sqrt(x1(2:)**2 + x2(2:)**2) * cos(2*sll_pi*dt + atan2(x2(2:),x1(2:)))
        x2_char(2:) = sqrt(x1(2:)**2 + x2(2:)**2) * sin(2*sll_pi*dt + atan2(x2(2:),x1(2:)))
     end if

     ! Time loop
     nloops = 0


     call cpu_time(t_init)
     print*,""
     do while (t .lt. tmax)
        !Error variables
        norm2_error = 0._f64
        diff_error  = 0._f64
        where_error = -1

        nloops = nloops + 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! let us compute the derivatives in every hexaedric direction with p the degree of the approximation
        ! then let us compute the derivatives in the x1 and x2 direction

        ! call der_finite_difference( f_tn, deriv, p, mesh ) -> CALCUL DE TOUTES LES DÉRIVÉES DANS LES DIRECTIONS H1, H2 et H3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        t      = t + dt

        do i=1, mesh%num_pts_tot   ! interpolation en chaque point

           ! ******************
           ! Approximation
           ! ******************

           ! computation of the interpolation   F( t_(n+1) , x_i , v_j ) = F ( t_n , X(t_n) , V(t_n) ) 

           !hermite_interpolation( x, y, f_tn, mesh, deriv, a, aire) 

           ! hermite_interpolation( f_tn , mesh, deriv, x1_char(i), x2_char(i) )  
           !-> contient toute l'interpolation à partir des 10 degrés de libertés


           ! if (f_tn(i) < 0.) then
           !    ! print *, "Negative value at"
           !    ! print *, "      Time  :", t
           !    ! print *, "      Loop  :", nloops
           !    ! print *, "      Point :", i
           !    ! print *, "      X1char:", x1_char(i)
           !    ! print *, "      X2char:", x2_char(i)
           !    ! STOP
           !    f_tn(i) = 0._f64
           ! end if

           ! ******************
           ! Analytical value    (-> in order to compute the error )
           ! ******************

           ! Computing characteristics

           if (which_advec .eq. 0) then
              ! linear advection
              x1(i) = from_global_x1(mesh, i-1) - advec*dt*nloops
              x2(i) = from_global_x2(mesh, i-1) - advec*dt*nloops
           else
              ! Circular advection
              x1_temp = sqrt(x1(i)**2 + x2(i)**2) * cos(2*sll_pi*dt + atan2(x2(i),x1(i)))
              x2(i)   = sqrt(x1(i)**2 + x2(i)**2) * sin(2*sll_pi*dt + atan2(x2(i),x1(i)))
              x1(i)   = x1_temp
           end if

           f_fin(i) = exp(-0.5_f64*((x1(i)-gauss_x1)**2/gauss_sig**2 &
                + (x2(i)-gauss_x2)**2 / gauss_sig**2))
           if (exponent(f_fin(i)) .lt. -17) then
              f_fin(i) = 0._f64
           end if



           !x1_basis = change_basis_x1(mesh, spline, x1(i), x2(i))
           !x2_basis = change_basis_x2(mesh, spline, x1(i), x2(i))   -> subroutine ne marche plus

           !chi1(i) = chi_gen_val(x1_basis, x2_basis, 1)
           !chi2(i) = chi_gen_val(x1_basis, x2_basis, 2)   -> subroutine ne marche plus - > tte des splines -> faire la meme chose pour mes éléments finis


           if (get_hex_num(i) .lt. get_hex_num(mesh%num_pts_tot) - deg*2 - 1) then
              ! Relative error
              if (diff_error .lt. abs(f_fin(i) - f_tn(i)) ) then
                 diff_error = abs(f_fin(i) - f_tn(i))
                 where_error = i
              end if
              ! Norm2 error :
              norm2_error = norm2_error + abs(f_fin(i) - f_tn(i))**2
           end if
        end do

        ! Norm2 error :
        norm2_error = sqrt(norm2_error)

        call cpu_time(t_end)
        ! Printing error
        print*,"  nt =", nloops, "     | error_Linf = ", diff_error
        print*,"                       | at hex =", get_hex_num(where_error), where_error
        print*,"                       | error_L2   = ", norm2_error
        ! print*," Center error = ", f_fin(1)-f_tn(1)




        !WRITING ERROR REGARDING TIME STEP
        if (t .eq. dt) then
           open (unit=12,file="err_chi2_gauss_cstadv.txt", &
                action="write",status="replace")
           write (12, "(3(g13.3,1x))") t, diff_error, norm2_error
           close(12)
        else 
           open (unit=12,file="err_chi2_gauss_cstadv.txt", &
                action="write",status="old", position="append")
           write (12, "(3(g13.3,1x))") t, diff_error, norm2_error
           close(12)
        end if

        call int2string(nloops,filenum)
        filename2 = "./time_files/analytical/ana_dist"//trim(filenum)//".txt"
        filename  = "./time_files/numerical/num_dist"//trim(filenum)//".txt"
        print*,filename
        print*,filename2
        call write_field_hex_mesh(mesh, f_tn, filename)
        call write_field_hex_mesh(mesh, f_fin,filename2)

     end do


     call write_field_hex_mesh(mesh, f_tn, "final_dist.txt")
     call write_field_hex_mesh(mesh, chi1, "chi1.txt")
     call write_field_hex_mesh(mesh, chi2, "chi2.txt")
     call write_field_hex_mesh(mesh, f_fin, "an_dist.txt")
     print *,""
     print*," *    Final error  = ", diff_error, " *"


     ! !WRITING ERROR REGARDING NUMBER OF POINTS
     ! if (num_cells .eq. 10) then 
     !    !NEW FILE :
     !    open (unit=12,file="error_file.txt",action="write",&
     !         status="replace")
     !    write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
     !    close(12)
     ! else
     !    !WRITE
     !    open (unit=12,file="error_file.txt",action="write",&
     !         status="old", position="append") 
     !    write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
     !    close(12)
     ! end if


     ! !WRITING CPU TIME
     ! if (num_cells .eq. 10) then 
     !    !NEW FILE :
     !    open (unit=12,file="cpu_time.txt",action="write",&
     !         status="replace")
     !    write (12, "(3(G15.3,1x))") num_cells, mesh%num_pts_tot, t_end-t_init
     !    close(12)
     ! else
     !    !WRITE
     !    open (unit=12,file="cpu_time.txt",action="write",&
     !         status="old", position="append") 
     !    write (12, "(3(G15.3,1x))") num_cells, mesh%num_pts_tot, t_end-t_init
     !    close(12)
     ! end if


     SLL_DEALLOCATE_ARRAY(f_init,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn,ierr)
     SLL_DEALLOCATE_ARRAY(f_fin,ierr)
     SLL_DEALLOCATE_ARRAY(x1,ierr)
     SLL_DEALLOCATE_ARRAY(x2,ierr)
     SLL_DEALLOCATE_ARRAY(chi1,ierr)
     SLL_DEALLOCATE_ARRAY(chi2,ierr)
     SLL_DEALLOCATE_ARRAY(x1_char,ierr)
     SLL_DEALLOCATE_ARRAY(x2_char,ierr)
     SLL_DEALLOCATE(mesh,ierr)


  end do


  contains
    
    
    type :: deriv
       sll_real64 :: dh1, dh2, dh3                     
    end type deriv

    
    subroutine der_finite_difference( f_tn, deriv, p, mesh) !-> CALCUL DE TOUTES LES DÉRIVÉES DANS LES DIRECTIONS H1, H2 ET H3
      implicit none
      type(hex_mesh_2d), pointer           :: mesh
      type(derive)     , pointer           :: deriv
      sll_int32                            :: i, j, r, s, n, n1, n2, n3
      sll_real64, dimension(:),allocatable :: w, ff, fff  
      sll_real64, intent(in)               :: p   
      sll_real64                           :: h1 ,h2

      !****************************************************************************************
      !computation of the coefficients in fonction of the degree p of the approximation required
      !****************************************************************************************
      if ( (p/2)*2 == p ) then !if p is even
         r = -p/2       
         s = p/2
      else  !if p is odd
         r = -p/2 - 1 ! de-centered to the left      
         s = p/2 !+ 1 ! de-centered to the right
      endif

      n = p + 1

      SLL_ALLOCATE( w(r,s), ff(r,s), fff(r,s), ffff(r,s),ierr )

      compute_w_hermite(w,r,s) 

      !************************************************************************
      ! Computation of the derivatives in both hexaedric directions h1 and h2,
      ! then for both x and y directions for every points of the mesh
      !************************************************************************
      
      do i = 1, mesh%num_pts_tot 

         
         !find the required points in the h1, h2 and h3 directions 

         ! find the hexa coordinate from the numerotation

         !->  voir la fonction de laura find-machin (...)  ???

         h1 = 
         h2 = 

         do j = r,s
            
            ! faire attention que le point est toujours dans le mesh -> sinon on applique la condition de bords

            ! trouver le numéro n1 du point correspondant à ( h1 + j , h2     )
            ! trouver le numéro n2 du point correspondant à ( h1     , h2 + j )
            ! trouver le numéro n3 du point correspondant à ( h1 + j , h2 + j )

            ff(j)   = f_tn( n1 ) 

            fff(j)  = f_tn( n2 ) 

            ffff(j) = f_tn( n3 ) 

         enddo


         dh1 = 0._f64
         dh2 = 0._f64
         dh3 = 0._f64

         do j = r,s
            dh1 = dh1 + w(j) * ff(j)  
            dh2 = dh2 + w(j) * fff(j)
            dh3 = dh3 + w(j) * ffff(j)
         enddo
         
         deriv(i)%dfh1 = dh1
         deriv(i)%dfh2 = dh2
         deriv(i)%dfh3 = dh3 

      enddo


     SLL_DEALLOCATE_ARRAY deallocate( w, ff, fff, ffff )

    end subroutine der_finite_difference


    

    subroutine hermite_interpolation( x, y, f_tn, mesh, deriv, a, aire) 

      implicit none

      type(derivative)           :: deriv
      type(hex_mesh_2d)          :: mesh
      sll_real64,intent(in)      :: x, y, a, aire
      sll_real64                 :: x1, x2, x3, y1, y2, y3
      sll_real64                 :: phi, ksi1, ksi2, ksi3, l1, l2, l3
      sll_real64                 :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32
      sll_real64,dimension(:)    :: f_tn
      sll_real64,dimension(1:10) :: freedom, base
      sll_int32                  :: i


      ! find the triangle where (x,y) belong to

      call find_tri(x,y,i)

      ! get the ten degrees of freedom

      ! values at the vertexes of the triangle

      freedom(1) = f_tn( mesh(i)%point(1) )
      freedom(2) = f_tn( mesh(i)%point(2) )
      freedom(3) = f_tn( mesh(i)%point(3) )

      ! value at the center of the triangle

      freedom(4) = 1.0_f64 / 3.0_f64 * ( freedom(1) + freedom(2) + freedom(3) )

      ! values of the derivatives

      !dans le calcul précédents des dérivées on a arrangé l'ordre 

      !find_point()

      ! i1 is the indice for the point S1 of the triangle
      ! i2 is the indice for the point S2 of the triangle
      ! i3 is the indice for the point S3 of the triangle

      ! il faut calculer les dérivées suivant les sommets -> (h1,h2) ou (h1,h3) ou (h2,h3)

      ! cas (h1,h2) -> inutile en fin de compte
      ! dfx = ( deriv(i1)%dfh1 - deriv(i1)%dfh2 ) * sqrt(3_f64) * 0.5_f64
      ! dfy = ( deriv(i1)%dfh1 + deriv(i1)%dfh2 ) * 0.5_f64

      freedom(5)  = deriv(i1)%dfh1
      freedom(6)  = deriv(i1)%dfh2

      ! cas (h1,h3) 
      ! dfx = deriv(i2)%dfh1 * sqrt(3_f64) * 0.5_f64
      ! dfy = deriv(i2)%dfh1 * 0.5_f64 + deriv(i2)%dfh3 

      freedom(7)  = deriv(i2)%dfh1 
      freedom(8)  = deriv(i2)%dfh3

      ! cas (h2,h3) 
      ! dfx = deriv(i3)%dfh2 * (-sqrt(3_f64)) * 0.5_f64
      ! dfy = deriv(i3)%dfh2 * 0.5_f64 + deriv(i3)%dfh3 

      freedom(9)   = deriv(i3)%dfh2 
      freedom(10)  = deriv(i3)%dfh3

      l1   = 0.5_f64 * abs( (x2 - x)*(y3 - y) - (x3 - x)*(y2 - y) ) / aire 
      l2   = 0.5_f64 * abs( (x1 - x)*(y3 - y) - (x3 - x)*(y1 - y) ) / aire 
      l3   = 0.5_f64 * abs( (x1 - x)*(y2 - y) - (x2 - x)*(y1 - y) ) / aire 

      phi = l1*l2*l3

      ksi1 = l1**3 - phi
      ksi2 = l2**3 - phi
      ksi3 = l3**3 - phi

      ksi12 = l1**2*l2 + phi*0.5_f64
      ksi13 = l1**2*l3 + phi*0.5_f64
      ksi21 = l2**2*l1 + phi*0.5_f64
      ksi23 = l2**2*l3 + phi*0.5_f64
      ksi31 = l3**2*l1 + phi*0.5_f64
      ksi32 = l3**2*l2 + phi*0.5_f64

      base(1) = 3.0_f64 * l1**2 - 2.0_f64*ksi1 - 9.0_f64*phi
      base(2) = 3.0_f64 * l2**2 - 2.0_f64*ksi2 - 9.0_f64*phi
      base(3) = 3.0_f64 * l3**2 - 2.0_f64*ksi3 - 9.0_f64*phi
      base(4) = 27.0_f64 * phi
      base(5) = a * ( ksi12 - 1.5_f64*phi )
      base(6) = a * ( ksi13 - 1.5_f64*phi )
      base(7) = a * ( ksi21 - 1.5_f64*phi )
      base(8) = a * ( ksi23 - 1.5_f64*phi )
      base(9) = a * ( ksi31 - 1.5_f64*phi )
      base(10)= a * ( ksi32 - 1.5_f64*phi )

      f = 0.0_f64

      do i = 1,10
         f = f + freedom(i)*base(i)
      enddo


    end subroutine hermite_interpolation


end program test_hex_hermite














