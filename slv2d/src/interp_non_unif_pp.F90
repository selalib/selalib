module interp_non_unif_pp_class
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_interpolators.h"
  use used_precision
  use geometry_module
  use sll_constants
  !use clock
  implicit none
  private
  public :: new, interpole
  type, public :: interp_non_unif_pp
     type (geometry) :: geom
     integer  :: conservative_x_case,mesh_x_case,interp_x_case(2),mesh_interp_x_case
     integer  :: conservative_y_case,mesh_y_case,interp_y_case(2),mesh_interp_y_case
     integer  :: bloc_index_x(3),bloc_index_y(3)
     real(wp) :: bloc_coord_x(2),bloc_coord_y(2)
     real(wp) :: eps_mesh_x,eps_mesh_y
     integer  :: rho_x_case,rho_y_case
     real(wp),dimension(:),allocatable :: node_positions_x,node_positions_y
     real(wp),dimension(:),allocatable :: node_pos_x,node_pos_y
     
  end type interp_non_unif_pp
  interface new
     module procedure new_interp_non_unif_pp
  end interface
  interface interpole
     module procedure interpole_non_unif_pp,interpole_non_unif_ppdep
  end interface
contains
  subroutine new_interp_non_unif_pp(this,geom,iflag)
    type(interp_non_unif_pp), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag

    ! initialisation de geom
    this%geom = geom
    
    this%mesh_interp_x_case = 1   ! 1: piecewise linear 2: cubic splines based
    this%eps_mesh_x = 0._f64     ! relax cubic splines based mesh with uniform
    this%mesh_interp_y_case = 1   ! 1: piecewise linear 2: cubic splines based
    this%eps_mesh_y = 0._f64     ! relax cubic splines based mesh with uniform


    !parameters
    ! bloc_coord = (a,b) \subset (0,1) zone to refine
    ! bloc_index(1) = density of points for (0,a) = 1 for the moment
    ! bloc_index(2) = density of points for (a,b)
    ! bloc_index(3) = density of points for (b,1) = 1 for the moment
    this%bloc_coord_x(1) = 0.53125_f64
    this%bloc_coord_x(2) = 0.6875_f64!69_f64
    this%bloc_index_x(1) = 1
    this%bloc_index_x(2) = 2!32
    this%bloc_index_x(3) = 1
    
    this%bloc_coord_y(1) = 0.53125_f64
    this%bloc_coord_y(2) = 0.6875_f64!69_f64
    this%bloc_index_y(1) = 1
    this%bloc_index_y(2) = 2!32
    this%bloc_index_y(3) = 1

    this%conservative_x_case=0
    this%mesh_x_case=1
    this%interp_x_case=(/1,3/)   ! cubic splines on uniform mesh
    this%interp_x_case=(/2,17/)  ! Lagrange 17 on uniform mesh
    !interp_case=(/0,17/)  ! FD order  17 on uniform mesh
    !interp_case=(/0,6/)  ! FD order  6 on uniform mesh
    !interp_case=(/10,5/) ! FD order 17 on two grid
    !interp_case=(/10,6/) ! FD order 6 on two grid
    !interp_case=(/11,3/)  ! cubic splines on two grid
    this%interp_x_case=(/31,3/)  ! cubic non uniform splines

    this%conservative_y_case=0
    this%mesh_y_case=1
    this%interp_y_case=(/1,3/)   ! cubic splines on uniform mesh
    this%interp_y_case=(/2,17/)  ! Lagrange 17 on uniform mesh
    this%interp_y_case=(/31,3/)  ! cubic non uniform splines

    this%rho_x_case = 1 !1:trapezoidal formula 
    this%rho_y_case = 1 !1:trapezoidal formula 

    allocate(this%node_positions_x(this%geom%nx+1))
    allocate(this%node_positions_y(this%geom%ny+1))
    
    call initialize(this)


 end subroutine new_interp_non_unif_pp

  subroutine interpole_non_unif_pp(this,fin,fout,x,y) 
    type(interp_non_unif_pp), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille precedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destine a contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points dans les quels on veut 
    ! evaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout on trouve en sortie les valeurs de f(i,j) 
    ! dans les points x(i),y(i).
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr

    

  end subroutine interpole_non_unif_pp

  subroutine interpole_non_unif_ppdep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(interp_non_unif_pp), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx, depy 
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    logical  :: aff
    integer  :: i,Nc1,Nc2,bloc_index(3)
    real(wp) :: alpha,mean,bloc_coord(2)
    real(wp) :: eps_mesh
    real(wp),dimension(:),allocatable :: buf_f,node_positions
    
    integer  :: conservative_case,mesh_case,interp_case(2),mesh_interp_case
    integer  :: transpose_case
    integer  :: N_buf
    real(wp),dimension(:),allocatable :: buf

    ALLOCATE(node_positions(max(this%geom%nx,this%geom%ny)+1))
    
    do transpose_case=0,1
      if(transpose_case==0)then
        Nc1=this%geom%nx
        Nc2=this%geom%ny
        bloc_coord=this%bloc_coord_x
        bloc_index=this%bloc_index_x
        mesh_case=this%mesh_x_case
        interp_case=this%interp_x_case
        mesh_interp_case=this%mesh_interp_x_case
        conservative_case=this%conservative_x_case
        eps_mesh=this%eps_mesh_x
        alpha = depx/(real(this%geom%nx,f64)*this%geom%dx)
        node_positions(1:Nc1+1)=this%node_positions_x(1:Nc1+1)
      else
        Nc1=this%geom%ny
        Nc2=this%geom%nx
        bloc_coord=this%bloc_coord_y
        bloc_index=this%bloc_index_y
        mesh_case=this%mesh_y_case
        interp_case=this%interp_y_case
        mesh_interp_case=this%mesh_interp_y_case
        conservative_case=this%conservative_y_case
        eps_mesh=this%eps_mesh_y
        alpha = depy/(real(this%geom%ny,f64)*this%geom%dy)    
        node_positions(1:Nc1+1)=this%node_positions_y(1:Nc1+1)
      endif
      
      
      
!enforce use of uniform grid depending on interp_case 
!      if((interp_case(1)>=0).and.(interp_case(1)<10))then
!        bloc_index(1) = 1
!        bloc_index(2) = 1
!        bloc_index(3) = 1
!        mesh_interp_case = 1
!        eps_mesh = 0._f64  !or 1._f64
!      endif
!  
!  
!enforce use of two_grid depending on interp_case 
!      if((interp_case(1)>=10).and.(interp_case(1)<29))then
!        bloc_index(1) = 1
!        bloc_index(3) = 1
!        mesh_interp_case = 1
!        eps_mesh = 0._f64  
!      endif
      
	  ALLOCATE(buf_f(Nc1+1))
	  !ALLOCATE(node_positions(Nc1+1))    
!      call compute_bloc(bloc_coord,bloc_index,Nc1)
!      if(mesh_case==1)then
!        do i=1,Nc1+1
!          node_positions(i) = (real(i,f64)-1._f64)/real(Nc1,f64)
!        enddo
!      endif
!      if(mesh_case==2)then
!        call compute_mesh_from_bloc(bloc_coord,bloc_index,node_positions,&
!          &mesh_interp_case,eps_mesh)
!      endif
      if((interp_case(1)>=0).and.(interp_case(1)<10))then  
        call constant_advection_size(N_buf,interp_case,Nc1)  
        allocate(buf(0:N_buf-1))
        call constant_advection_init(buf,N_buf,Nc1,interp_case)
      endif
      do i=1,Nc2
        if(transpose_case==0)then
          buf_f(1:Nc1)=f(1:Nc1,i)
          buf_f(Nc1+1)=f(1,i)
        else
          buf_f(1:Nc1)=f(i,1:Nc1)
          buf_f(Nc1+1)=f(i,1)      
        endif  
        if(conservative_case==1)then
          call function_to_primitive(buf_f,node_positions,Nc1,mean)  
        endif       
        !choice of interpolation in v with different interfaces
        !uniform
        if((interp_case(1)>=0).and.(interp_case(1)<10))then
          call constant_advection_solve(buf,N_buf,buf_f(1:Nc1),Nc1,alpha,&
          &interp_case)
          buf_f(Nc1+1)=buf_f(1)
        endif
        !advective two_grid
        if((interp_case(1)>=10).and.(interp_case(1)<20))then
          call constant_advection_two_grid_per(buf_f,alpha,&
          &bloc_coord,bloc_index,node_positions,interp_case)
        endif
        !advective non_uniform
        if((interp_case(1)>=30).and.(interp_case(1)<40))then
          call constant_advection_non_unif_per(buf_f,alpha,node_positions,&
          &Nc1,interp_case)
        endif
        if(conservative_case==1)then
          call primitive_to_function(buf_f,node_positions,Nc1,mean)  
        endif
        if(transpose_case==0)then
          f(1:Nc1,i)=buf_f(1:Nc1)
        else
          f(i,1:Nc1)=buf_f(1:Nc1)        
        endif  
      enddo
      if((interp_case(1)>=0).and.(interp_case(1)<10))then  
        DEALLOCATE(buf)
      endif
      DEALLOCATE(buf_f)
    enddo
    !if (aff) then 
    !   call clck_temps(l_a)
    !end if
    DEALLOCATE(node_positions)
      
  end subroutine interpole_non_unif_ppdep



  subroutine initialize(this) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(interp_non_unif_pp), intent(inout) :: this
    integer  :: i,Nc1,Nc2,bloc_index(3)
    real(wp) :: alpha,mean,bloc_coord(2)
    real(wp) :: eps_mesh
    real(wp),dimension(:),allocatable :: buf_f,node_positions
    
    integer  :: conservative_case,mesh_case,interp_case(2),mesh_interp_case
    integer  :: transpose_case
    integer  :: N_buf
    real(wp),dimension(:),allocatable :: buf

    
    do transpose_case=0,1
      if(transpose_case==0)then
        Nc1=this%geom%nx
        Nc2=this%geom%ny
        bloc_coord=this%bloc_coord_x
        bloc_index=this%bloc_index_x
        mesh_case=this%mesh_x_case
        interp_case=this%interp_x_case
        mesh_interp_case=this%mesh_interp_x_case
        conservative_case=this%conservative_x_case
        eps_mesh=this%eps_mesh_x
      else
        Nc1=this%geom%ny
        Nc2=this%geom%nx
        bloc_coord=this%bloc_coord_y
        bloc_index=this%bloc_index_y
        mesh_case=this%mesh_y_case
        interp_case=this%interp_y_case
        mesh_interp_case=this%mesh_interp_y_case
        conservative_case=this%conservative_y_case
        eps_mesh=this%eps_mesh_y
      endif
      
      
      
      !enforce use of uniform grid depending on interp_case 
      if((interp_case(1)>=0).and.(interp_case(1)<10))then
        bloc_index(1) = 1
        bloc_index(2) = 1
        bloc_index(3) = 1
        mesh_interp_case = 1
        eps_mesh = 0._f64  !or 1._f64
      endif
  
  
      !enforce use of two_grid depending on interp_case 
      if((interp_case(1)>=10).and.(interp_case(1)<29))then
        bloc_index(1) = 1
        bloc_index(3) = 1
        mesh_interp_case = 1
        eps_mesh = 0._f64  
      endif
      
      ALLOCATE(node_positions(Nc1+1))
      if(transpose_case==0)then
        node_positions=this%node_positions_x
      else
        node_positions=this%node_positions_y      
      endif
      call compute_bloc(bloc_coord,bloc_index,Nc1)
      if(mesh_case==1)then
        do i=1,Nc1+1
          node_positions(i) = (real(i,f64)-1._f64)/real(Nc1,f64)
        enddo
      endif
      if(mesh_case==2)then
        call compute_mesh_from_bloc(bloc_coord,bloc_index,node_positions,&
          &mesh_interp_case,eps_mesh)
      endif
      
      if(transpose_case==0)then
        this%node_positions_x=node_positions
      else
        this%node_positions_y=node_positions     
      endif
      
      !initialize node_pos
      if(conservative_case==1)then 
        do i=1,Nc1
          node_positions(i)= 0.5_f64*(node_positions(i)+node_positions(i+1))
        enddo
        node_positions(Nc1+1)=node_positions(1)+1._f64
      endif

      if(transpose_case==0)then
        this%node_pos_x=node_positions
      else
        this%node_pos_y=node_positions     
      endif
      
      
      DEALLOCATE(node_positions)
    enddo
    !if (aff) then 
    !   call clck_temps(l_a)
    !end if
      
  end subroutine initialize





!!!!!!!!!!!!!! routines from vlasov1d_non_unif_v_module.F90

  
  subroutine compute_dual_mesh(node_positions,node_positions_dual,N)
    sll_real64,dimension(:),allocatable::node_positions,node_positions_dual
    sll_int32,intent(in)::N
    sll_int32::i
    do i=2,N+1
      node_positions_dual(i)=0.5_f64*(node_positions(i)+node_positions(i-1))
    enddo
    node_positions_dual(1)=0.5_f64*(node_positions(1)+node_positions(N)-1._f64)
    node_positions_dual=node_positions_dual-node_positions_dual(1)

  end subroutine compute_dual_mesh
  
  subroutine compute_bloc(bloc_coord,bloc_index,N)
    !input: a=bloc_coord(1) b=bloc_coord(2)
    !       (a,b) \subset (0,1) is the refine zone
    !       bloc_index(1) = density of points in (0,a) 
    !       bloc_index(2) = density of points in (a,b) 
    !       bloc_index(3) = density of points in (b,1)
    !output:0<=i1<i1+N_fine<=N and x(i1)=a, x(i1+N_fine)=b (approx), x(0)=0, x(N)=1
    !       bloc_coord(1) = x(i1)
    !       bloc_coord(2) = x(i1+N_fine)
    !       bloc_index(1) = i1 
    !       bloc_index(2) = N_fine 
    !       bloc_index(3) = N-i1-N_fine
    sll_real64,intent(inout)::bloc_coord(2)
    sll_int32,intent(inout)::bloc_index(3)
    sll_int32,intent(in)::N
    sll_real64::a,b
    sll_int32::i1,i2,N_coarse,N_local,N_fine
    
    a=bloc_coord(1)
    b=bloc_coord(2)
    
    
    !case of uniform mesh with refined zone
    !we have a coarse mesh with N_coarse
    !N=i1+N_local*(i2-i1)+N_coarse-i2
    !N_fine=N_local*(i2-i1)
    !x(i1)=i1/N_coarse x(i1+N_fine)=i2/N_coarse
    if((bloc_index(1)==1).and.(bloc_index(3)==1))then      
      N_local = bloc_index(2)
      N_coarse = floor(real(N,f64)/(1._f64+(b-a)*(real(N_local,f64)-1._f64)))
      if(N_local/=1)then
        i2 = (N-N_coarse)/(N_local-1)
      else
        i2 = floor((b-a)*N_coarse)  
      endif   
      N_coarse = N-i2*(N_local-1)
      i1 = floor(a*N_coarse)
      i2 = i2+i1
      bloc_index(1)=i1
      N_fine=N_local*(i2-i1)
      bloc_index(2)=N_fine
      bloc_index(3)=N-i1-N_fine
      bloc_coord(1)=real(i1,f64)/real(N_coarse,f64)
      bloc_coord(2)=real(i2,f64)/real(N_coarse,f64)
           
      !print *,'#uniform fine mesh would be:',N_coarse*N_local
      !print *,'#N_coarse=',N_coarse
      !print *,'#saving:',real(N,f64)/real(N_coarse*N_local,f64)
      !print *,'#new x(i1),x(i1+N_fine)=',bloc_coord(1),bloc_coord(2)
      !print *,'#error for x(i1),x(i1+N_fine)=',bloc_coord(1)-a,bloc_coord(2)-b
      !print *,'#i1,i1+N_fine,N_fine,N=',i1,i1+N_fine,N_fine,N
    else
      print *,'case in compute_bloc not implemented yet'
      stop  
    endif
    
    
  end subroutine compute_bloc
  
  
  subroutine compute_mesh_from_bloc(bloc_coord,bloc_index,node_positions,&
  &mesh_interp_case,eps_mesh)
    !input:   x1=bloc_coord(1),x2=bloc_coord(2)
    !         with 0<i1<i2<N
    !         i1=bloc_index(1), i2=i1+bloc_index(2)
    !         N=bloc_index(1)+bloc_index(2)+bloc_index(3)
    !         mesh_interp_case !1: linear 2: cubic splines
    !         eps_mesh ! uniform relation for cubic splines mesh
    !output:  node_positions(1:N+1) depending on mesh_interp_case
    !         with constraints node_positions(i1+1)=x1,node_positions(i2+1)=x2
    !         node_positions(1)=0, node_positions(N+1)=1
    sll_int32,intent(in)::bloc_index(3),mesh_interp_case
    sll_real64,intent(in)::bloc_coord(2),eps_mesh
    sll_real64,dimension(:),allocatable::node_positions
    sll_int32::i,i1,i2,N
    sll_real64::dx
    
    N=bloc_index(1)+bloc_index(2)+bloc_index(3)
    i1=bloc_index(1)
    i2=i1+bloc_index(2)
    node_positions(1:N+1)=-1._f64
    node_positions(1)=0._f64
    node_positions(i1+1)=bloc_coord(1)
    node_positions(i2+1)=bloc_coord(2)
    node_positions(N+1)=1._f64
    
    
    !piecewise linear mapping (maybe enhanced like in complete mesh)
    if(mesh_interp_case==1)then
      dx=bloc_coord(1)/real(bloc_index(1),f64)
      do i=2,bloc_index(1)
        node_positions(i) = (real(i,f64)-1._f64)*dx
      enddo
      dx=(bloc_coord(2)-bloc_coord(1))/real(bloc_index(2),f64)
      do i=2,bloc_index(2)
        node_positions(i+i1)=bloc_coord(1)+(real(i,f64)-1._f64)*dx
      enddo
      dx=(1._f64-bloc_coord(2))/real(bloc_index(3),f64)
      do i=2,bloc_index(3)
        node_positions(i+i2)=bloc_coord(2)+(real(i,f64)-1._f64)*dx
      enddo
    endif
    if(mesh_interp_case==2)then
      call complete_mesh_unit(node_positions,N)      
      do i=2,N
        node_positions(i)=&
          &(1-eps_mesh)*node_positions(i)+eps_mesh*(real(i,f64)-1._f64)/real(N,f64)
      enddo
    endif
    
    if((mesh_interp_case/=1).and.(mesh_interp_case/=2))then
      print *,'#bad mesh_interp_case=',mesh_interp_case
      print *,'#in subroutine compute_mesh_from_bloc'
      stop
    endif
          
  end subroutine compute_mesh_from_bloc
  
  
  
  subroutine complete_mesh_unit(x,N)
    sll_real64,dimension(:),allocatable::x
    sll_int32,intent(in)::N
    !we have the mesh x(1:N+1), with x(1)=0, x(N+1)=1
    !if k is such that x(k)=-1, we have to put a new value
    !obtained here with cubic non uniform splines
    sll_int32::i,N_cell,ii
    sll_real64::tmp
    !temporary allocations
    sll_real64,dimension(:),pointer :: buf,Xstar,node_pos,coeffs
    sll_real64,dimension(:),allocatable::f
    sll_int32,dimension(:),pointer :: ibuf 

    
    !we first count the numbers of cells
    !and check that the input is valid
    if((x(1)/=0._f64).or.(x(N+1)/=1._f64))then
      print *,'bad value for (x(1),x(N+1))=',x(1),x(N+1)
      print *,'should be (0,1)'
      stop
    endif      
    N_cell=0
    tmp=x(1)
    do i=2,N+1
      if(x(i)/=-1._f64)then
        if(x(i)<=tmp)then
          print *,'x(',i,')=',x(i),' should be greater than',tmp
          stop
        endif
        N_cell=N_cell+1
        tmp=x(i)
      endif
    enddo
    
    print *,'#N_cell=',N_cell

    allocate(buf(10*N_cell))
    allocate(ibuf(N_cell))
    allocate(node_pos(-1:N_cell+2),coeffs(-1:N_cell+2))
    allocate(Xstar(1:N+1))
    allocate(f(1:N_cell+1))
    
    !fill node_pos array 
    ii=0
    do i=1,N+1
      if(x(i)/=-1._f64)then
        node_pos(ii)=(real(i,f64)-1._f64)/real(N,f64)
        f(ii+1)=x(i)-node_pos(ii)   
        ii=ii+1
      endif
    enddo
    
    do i=1,N+1
      Xstar(i)=(real(i,f64)-1._f64)/real(N,f64)
    enddo
    
    !do ii=0,N_cell
    !  print *,ii,node_pos(ii),f(ii+1)
    !enddo
    
        
    call setup_spline_nonunif_1D_periodic_aux( node_pos, N_cell, buf, ibuf)
    call compute_spline_nonunif_1D_periodic_aux( f, N_cell, buf, ibuf, coeffs )
    call interpolate_array_value_nonunif_aux&
         &( Xstar(1:N+1), x(1:N+1), N+1, node_pos, coeffs,N_cell)
    
    do i=1,N+1
      x(i)=x(i)+Xstar(i)
    enddo
    x(1)=0._f64
    x(N+1)=1._f64
    
    !do i=1,N+1
    !  print *,x(i),Xstar(i)
    !enddo
    
    !stop
  
  end subroutine complete_mesh_unit
  
  subroutine constant_advection_hermite_non_unif_per(f,df,alpha,node_positions,N)
    !alpha and node_positions are normalized to [0,1]
    implicit none
    sll_real64,intent(in)::alpha
    sll_int32,intent(in)::N
    sll_real64,dimension(:),allocatable::f
    sll_real64,dimension(:,:),allocatable::df
    sll_real64,dimension(:),allocatable::node_positions
    sll_real64,dimension(:),pointer :: Xstar
    sll_real64,dimension(:),allocatable :: f_old
    sll_real64 :: x,x_old,xx,w(0:1,0:1)
    sll_int32  :: i,j
        
    allocate(Xstar(1:N+1))
    allocate(f_old(1:N+1))
    
    do i=1,N+1    
      Xstar(i) = node_positions(i)-alpha
    enddo
    do i=1,N+1
      do while (Xstar(i).gt.1._f64)
        Xstar(i) = Xstar(i)-1._f64
      end do
      do while (Xstar(i).lt.0._f64)
        Xstar(i) = Xstar(i)+1._f64
      end do    
    enddo
    
    !do i=1,N+1
    !  print *,i,Xstar(i),node_positions(i)
    !enddo
    
    
    f_old = f
    x_old = Xstar(1)
    j=1
    do i=1,N+1 
      x = Xstar(i)
      if (x==1.0_f64) then
        j = N
        xx = 1._f64  
      else
        if(x.ge.x_old) then
          do while(node_positions(j+1).le.x)
            j = j+1
          enddo
        else
          do while(node_positions(j+1).gt.x)
            j = j-1
          enddo
          j=j+1
        endif
        !print *,i,j
        !stop
        if( .not.( (x.ge.node_positions(j)).and.(x.lt.node_positions(j+1)) ) ) then
          print *,'bad localization in constant_advection_hermite_non_unif_per'
          print *,'j=',j,'N=',N        
          print *,node_positions(j),'<=',x,'<',node_positions(j+1)
          stop
        endif
        xx=(x-node_positions(j))/(node_positions(j+1)-node_positions(j))
      endif
      x_old=x
      w(1,0)=xx**2*(3._f64-2._f64*xx)
      w(0,0)=1._f64-w(1,0)
      w(0,1)=xx*(1._f64-xx)**2
      w(1,1)=xx**2*(xx-1._f64)
      f(i)=w(0,0)*f_old(j)+w(1,0)*f_old(j+1)+w(0,1)*df(1,j)+w(1,1)*df(2,j)
      !f(i)=(1._f64-xx)*f_old(j)+xx*f_old(j+1)
        
      
    enddo
    
    
    deallocate(Xstar)
    deallocate(f_old)
    
  end subroutine constant_advection_hermite_non_unif_per


  subroutine constant_advection_two_grid_per(f,alpha,bloc_coord,bloc_index,node_positions,interp_case)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    sll_real64,dimension(:,:),allocatable::df
    sll_int32,intent(in):: bloc_index(3),interp_case(2)
    sll_real64,intent(in)::alpha,bloc_coord(2)
    sll_int32::i,N
    
    N=bloc_index(1)+bloc_index(2)+bloc_index(3)
    
    if(interp_case(1)==10)then
      allocate(df(2,N+1))
      df=0._f64
      call compute_hermite_derivatives_two_grid(f,df,bloc_coord,bloc_index,interp_case(2))
      call constant_advection_hermite_non_unif_per(f,df,alpha,node_positions,N)
      deallocate(df)
    endif    
    if(interp_case(1)==11)then
      allocate(df(2,N+1))
      df=0._f64
      if(interp_case(2)==3)then
        call compute_hermite_derivatives_two_grid_spl(f,df,bloc_coord,bloc_index)
      endif      
      call constant_advection_hermite_non_unif_per(f,df,alpha,node_positions,N)      
      deallocate(df)
    endif    


    !if(interp_case(1)==20)then
    !  call constant_advection_spl_non_unif_per(f,alpha,node_positions,N)
    !endif
    
    if((interp_case(1)/=10).and.((interp_case(1)/=11)))then
      print*,'#interp_case(1)=',interp_case(1),'not implemented'
      print *,'#in subroutine constant_advection_two_grid_per'
    endif
    
  end subroutine constant_advection_two_grid_per




  
  subroutine constant_advection_non_unif_per(f,alpha,node_positions,N,interp_case)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    !sll_int32,intent(in)::igeom(4)
    !sll_real64,dimension(:,:),allocatable::df
    sll_int32,intent(in):: N,interp_case(2)
    sll_real64,intent(in)::alpha
    
!    if(interp_case(1)==10)then
!      allocate(df(2,N+1))
!      df=0._f64
!      call compute_hermite_derivatives_piecewise_uniform(f,df,igeom,N,interp_case(2))
!      call constant_advection_hermite_non_unif_per(f,df,alpha,node_positions,N)
!      deallocate(df)
!    endif    
    if(interp_case(1)==31)then
      if(interp_case(2)==3)then
        call constant_advection_spl_non_unif_per(f,alpha,node_positions,N)
      endif  
    endif
    
    if((interp_case(1)/=31).and.((interp_case(1)/=31)))then
      print*,'#interp_case(1)=',interp_case(1),'not implemented'
      print *,'#in subroutine constant_advection_non_unif_per'
    endif
    
  end subroutine constant_advection_non_unif_per




  subroutine csl_constant_advection_non_unif_per(f,alpha,node_positions,N,interp_case)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    sll_int32,intent(in):: N,interp_case(2)
    sll_real64,intent(in)::alpha
    
    
    
    
    if(interp_case(1)==41)then
      if(interp_case(2)==3)then
        call csl_constant_advection_spl_non_unif_per(f,alpha,node_positions,N)
      endif  
    endif
    if((interp_case(1)/=41).and.((interp_case(1)/=41)))then
      print*,'#interp_case(1)=',interp_case(1),'not implemented'
      print *,'#in subroutine csl_constant_advection_non_unif_per'
    endif
    
  end subroutine csl_constant_advection_non_unif_per

  subroutine compute_hermite_derivatives_two_grid_spl(f,df,bloc_coord,bloc_index)
    implicit none
    sll_real64,dimension(:),allocatable::f
    sll_real64,intent(in)::bloc_coord(2)
    sll_real64,dimension(:,:),allocatable::df
    sll_int32,intent(in):: bloc_index(3)
    sll_int32::r,s,i,j,ii,N,N_coarse,N_local,rd,qd,jj,igeom(4)
    sll_real64::tmp,f_tmp(2,2),ww(0:1,0:1),xx
    sll_real64,dimension(:),allocatable::buf,buf2
    
    N=bloc_index(1)+bloc_index(2)+bloc_index(3)
    igeom(1)=bloc_index(1)
    tmp=real(bloc_index(1),f64)/bloc_coord(1)
    N_coarse = floor(tmp)
    !print *,tmp-N_coarse,tmp-N_coarse-1
    if(abs(tmp-N_coarse-1)<abs(tmp-N_coarse))then
      N_coarse=N_coarse+1
    endif
    !print *,N_coarse
    !print *,real(bloc_index(1),f64)/bloc_coord(1)
    !print *,'N_coarse (in compute_hermite_derivatives_two_grid_spl)=',N_coarse
    !stop
    igeom(2)=N_coarse-bloc_index(3)
    N_local=bloc_index(2)/(igeom(2)-igeom(1))
    
    igeom(3)=N_coarse
    igeom(4)=N_local
    
    
    allocate(buf(1:N_coarse+1))
    allocate(buf2(-1:N_local*(igeom(2)-igeom(1))+1))
    !we first compute the derivatives on the first unrefined zone
    do i=1,igeom(1)+1
      buf(i)=f(i)  
    enddo
    ii=igeom(1)+1
    do i=igeom(1)+2,igeom(2)+1     
      ii=ii+N_local
      buf(i)=f(ii)
    enddo
    do i=igeom(2)+2,N_coarse+1
      ii=ii+1
      buf(i)=f(ii)
    enddo
    !i1=igeom(1) i2=igeom(2)
    ! (1..i1+1) (i1+1..i2+1) (i2+1..N_coarse+1)
    
    !here we have to compute the derivatives with periodic splines
    ! to complete inout is buf(1:N+1)
    call compute_derivative_spl_per(buf,N_coarse)
    !!we compute on the refined middle zone
    buf2(-1)=buf(igeom(1)+1)/real(N_local,f64)
    do i=0,N_local*(igeom(2)-igeom(1))
      buf2(i)=f(igeom(1)+1+i)    
    enddo
    buf2(N_local*(igeom(2)-igeom(1))+1)=buf(igeom(2)+1)/real(N_local,f64)
    !compute the derivatives with hermite splines  
    !to complete inout is  buf2(-1:N_local*(igeom(2)-igeom(1))+1)
    !result in buf2(1:N_local*(igeom(2)-igeom(1))-1)
    call compute_derivative_spl_hrmt(buf2,N_local*(igeom(2)-igeom(1)))
    !copy to df
    do i=1,igeom(1)+1
      df(1,i)=buf(i)  
    enddo
    do i=1,N_local*(igeom(2)-igeom(1))-1 !interior points      
    !do i=igeom(1)+2,igeom(1)+2+N_local*(igeom(2)-igeom(1))!igeom(2)+1     
      df(1,i+igeom(1)+1)=buf2(i)
    enddo
    j=igeom(1)+N_local*(igeom(2)-igeom(1))
    do i=igeom(2)+1,N_coarse+1
      j=j+1
      df(1,j)=buf(i) 
    enddo
    do i=1,N
      df(2,i)=df(1,i+1)
    enddo
    df(2,N+1)=df(2,1)          
    df(1,igeom(1)+1)=df(2,igeom(1))/real(N_local,f64)
    j=igeom(1)+N_local*(igeom(2)-igeom(1))
    df(2,j)=df(1,j+1)/real(N_local,f64)

    deallocate(buf,buf2)
    
    !df=0._f64
     
  end subroutine compute_hermite_derivatives_two_grid_spl

  
  subroutine compute_hermite_derivatives_two_grid(f,df,bloc_coord,bloc_index,p)
    implicit none
    sll_real64,dimension(:),allocatable::f
    sll_real64,intent(in)::bloc_coord(2)
    sll_real64,dimension(:,:),allocatable::df
    sll_int32,intent(in):: bloc_index(3),p    
    sll_int32::r,s,i,j,ii,N,N_coarse,N_local,rd,qd,jj,igeom(4)
    sll_real64::tmp,f_tmp(2,2),ww(0:1,0:1),xx
    sll_real64,dimension(:),allocatable::buf,buf2
    sll_real64,dimension(:),allocatable::w
    
    N=bloc_index(1)+bloc_index(2)+bloc_index(3)
    !igeom(1)=i1, igeom(2)=i2
    
    igeom(1)=bloc_index(1)
    tmp=real(bloc_index(1),f64)/bloc_coord(1)
    N_coarse = floor(tmp)
    !print *,tmp-N_coarse,tmp-N_coarse-1
    if(abs(tmp-N_coarse-1)<abs(tmp-N_coarse))then
      N_coarse=N_coarse+1
    endif
    igeom(2)=N_coarse-bloc_index(3)
    N_local=bloc_index(2)/(igeom(2)-igeom(1))
    
    igeom(3)=N_coarse
    igeom(4)=N_local
    
    r=-p/2
    s=(p+1)/2
    allocate(w(r:s))
    
    
    
    allocate(buf(1:N_coarse+1))
    allocate(buf2(r:N_local*(igeom(2)-igeom(1))+s))
    
    call compute_finite_difference_init(w,p)
    

    !we first compute the derivatives on the first unrefined zone
    do i=1,igeom(1)+1
      buf(i)=f(i)  
    enddo
    ii=igeom(1)+1
    do i=igeom(1)+2,igeom(2)+1     
      ii=ii+N_local
      buf(i)=f(ii)
    enddo
    do i=igeom(2)+2,N_coarse+1
      ii=ii+1
      buf(i)=f(ii)
    enddo    
    
    
    !do i=1,N_coarse+1
    !  print *,i,buf(i)
    !enddo
    !stop
    
    do i=1,igeom(1)
      tmp=0._f64
      do ii=r,s
        tmp=tmp+w(ii)*buf(modulo(i+ii-1+N_coarse,N_coarse)+1)
      enddo
      df(1,i)=tmp
      tmp=0._f64
      if(-r/=s)then
        do ii=r,s
          tmp=tmp-w(r+s-ii)*buf(modulo(i+ii-1+N_coarse,N_coarse)+1)
        enddo
      else
        do ii=r,s
          tmp=tmp-w(r+s-ii)*buf(modulo(i+ii+N_coarse,N_coarse)+1)
        enddo        
      endif  
      df(2,i)=tmp
    enddo
    
    
    
    !we compute now on the last third unrefined zone
    j=igeom(1)+N_local*(igeom(2)-igeom(1))
    do i=igeom(2)+1,N_coarse
      j=j+1
      tmp=0._f64
      do ii=r,s
        tmp=tmp+w(ii)*buf(modulo(i+ii-1+N_coarse,N_coarse)+1)
      enddo
      df(1,j)=tmp
      tmp=0._f64
      if(-r/=s)then
        do ii=r,s
          tmp=tmp-w(r+s-ii)*buf(modulo(i+ii-1+N_coarse,N_coarse)+1)
        enddo
      else
        do ii=r,s
          tmp=tmp-w(r+s-ii)*buf(modulo(i+ii+N_coarse,N_coarse)+1)
        enddo      
      endif  
      df(2,j)=tmp
    enddo
    
    !!we compute on the refined middle zone
    do i=0,N_local*(igeom(2)-igeom(1))
      buf2(i)=f(igeom(1)+1+i)    
    enddo
    !we have to compute buf2(-1),..buf2(r)
    !(-r)=N_local*qd+rd, 0<=rd<qd
    qd=(-r)/N_local
    rd=(-r)-N_local*qd    
    !print *,qd,rd    
    i=0
    do j=1,qd
      f_tmp(1,1)=f(igeom(1)+1-j)!f_loc(0)
      f_tmp(2,1)=f(igeom(1)+2-j)!f_loc(1)
      f_tmp(1,2)=df(1,igeom(1)+1-j)!df_loc(0+)
      f_tmp(2,2)=df(2,igeom(1)+1-j)!df_loc(1-)      
      do jj=1,N_local
        i=i-1
        xx=1._f64-real(jj,f64)/real(N_local,f64)
        ww(1,0)=xx**2*(3._f64-2._f64*xx)
        ww(0,0)=1._f64-ww(1,0)
        ww(0,1)=xx*(1._f64-xx)**2
        ww(1,1)=xx**2*(xx-1._f64)
        buf2(i)=ww(0,0)*f_tmp(1,1)+ww(1,0)*f_tmp(2,1)+ww(0,1)*f_tmp(1,2)+ww(1,1)*f_tmp(2,2)
      enddo  
    enddo
    j=qd+1
    f_tmp(1,1)=f(igeom(1)+1-j)!f_loc(0)
    f_tmp(2,1)=f(igeom(1)+2-j)!f_loc(1)
    f_tmp(1,2)=df(1,igeom(1)+1-j)!df_loc(0+)
    f_tmp(2,2)=df(2,igeom(1)+1-j)!df_loc(1-)      
    do jj=1,rd
      i=i-1
      xx=1._f64-real(jj,f64)/real(N_local,f64)
      ww(1,0)=xx**2*(3._f64-2._f64*xx)
      ww(0,0)=1._f64-ww(1,0)
      ww(0,1)=xx*(1._f64-xx)**2
      ww(1,1)=xx**2*(xx-1._f64)
      buf2(i)=ww(0,0)*f_tmp(1,1)+ww(1,0)*f_tmp(2,1)+ww(0,1)*f_tmp(1,2)+ww(1,1)*f_tmp(2,2)      
    enddo    
    !print*,i,r
    i=N_local*(igeom(2)-igeom(1))
    !we have to compute buf2(i+1),..buf2(i+s-1)
    !s-1=N_local*qd+rd, 0<=rd<qd
    !qd=(s-1)/N_local
    !rd=(s-1)-N_local*qd    

    qd=(s)/N_local
    rd=(s)-N_local*qd    


    !to change
    ii=igeom(1)+N_local*(igeom(2)-igeom(1))
    do j=1,qd
      f_tmp(1,1)=f(ii+j)!f_loc(0)
      f_tmp(2,1)=f(ii+j+1)!f_loc(1)
      f_tmp(1,2)=df(1,ii+j)!df_loc(0+)
      f_tmp(2,2)=df(2,ii+j)!df_loc(1-)      
      do jj=1,N_local
        i=i+1
        xx=real(jj,f64)/real(N_local,f64)
        ww(1,0)=xx**2*(3._f64-2._f64*xx)
        ww(0,0)=1._f64-ww(1,0)
        ww(0,1)=xx*(1._f64-xx)**2
        ww(1,1)=xx**2*(xx-1._f64)
        buf2(i)=ww(0,0)*f_tmp(1,1)+ww(1,0)*f_tmp(2,1)+ww(0,1)*f_tmp(1,2)+ww(1,1)*f_tmp(2,2)
      enddo  
    enddo
    j=qd+1
    f_tmp(1,1)=f(ii+j)!f_loc(0)
    f_tmp(2,1)=f(ii+j+1)!f_loc(1)
    f_tmp(1,2)=df(1,ii+j)!df_loc(0+)
    f_tmp(2,2)=df(2,ii+j)!df_loc(1-)      
    do jj=1,rd
      i=i+1
      xx=real(jj,f64)/real(N_local,f64)
      ww(1,0)=xx**2*(3._f64-2._f64*xx)
      ww(0,0)=1._f64-ww(1,0)
      ww(0,1)=xx*(1._f64-xx)**2
      ww(1,1)=xx**2*(xx-1._f64)
      buf2(i)=ww(0,0)*f_tmp(1,1)+ww(1,0)*f_tmp(2,1)+ww(0,1)*f_tmp(1,2)+ww(1,1)*f_tmp(2,2)      
    enddo    
    !print*,i,N_local*(igeom(2)-igeom(1))+s-1
    
    !do i=r,N_local*(igeom(2)-igeom(1))+s-1
    !  print *,i,buf2(i)
    !enddo
    !stop
    !now, we have buf2 and we can compute
    !the derivatives in the middle refined zone


    j=igeom(1)
    do i=0,N_local*(igeom(2)-igeom(1))-1
      j=j+1
      tmp=0._f64
      do ii=r,s
        tmp=tmp+w(ii)*buf2(i+ii)
      enddo
      df(1,j)=tmp
      tmp=0._f64
      if(-r/=s)then
        do ii=r,s
          tmp=tmp-w(r+s-ii)*buf2(i+ii)
        enddo
      else
        do ii=r,s
          tmp=tmp-w(r+s-ii)*buf2(i+1+ii)
        enddo
      endif  
      df(2,j)=tmp
    enddo
    
    !do i=1,N+1
    !  print *,i,df(1,i),df(2,i)
    !enddo

    do i=1,igeom(1)
      !print *,(real(i,f64)-1._f64)/real(N_coarse,f64),df(1,i)*real(N_coarse,f64),df(2,i)*real(N_coarse,f64)
    enddo
    
    ii=igeom(1)*N_local+1
    do i=igeom(1)+1,igeom(1)+(igeom(2)-igeom(1))*N_local
      ii=ii+1
      !print *,(real(ii,f64)-1._f64)/real(N_coarse*N_local,f64),&
      !  &df(1,i)*real(N_coarse*N_local,f64),&
      !  &df(2,i)*real(N_coarse*N_local,f64)
    enddo
    
    ii=igeom(2) 
    do i=igeom(1)+(igeom(2)-igeom(1))*N_local+1,N
      ii=ii+1
      !print *,(real(ii,f64)-1._f64)/real(N_coarse,f64),df(1,i)*real(N_coarse,f64),df(2,i)*real(N_coarse,f64)
    enddo



    !stop
    deallocate(buf,buf2,w)
    
    
  end subroutine compute_hermite_derivatives_two_grid
  
  
  
  subroutine constant_advection_spl_non_unif_per(f,alpha,node_positions,N)
    !alpha and node_positions are normalized to [0,1]
    !use numeric_constants
    !use cubic_non_uniform_splines
    implicit none
    
    sll_real64,dimension(:),allocatable::f,node_positions
    !type(cubic_nonunif_spline_1D), pointer :: spl_per
    sll_int32,intent(in):: N
    sll_real64,intent(in)::alpha
    sll_real64 :: dx
    sll_int32  :: i
    sll_real64 :: M,tmp,tmp2
    !temporary allocations
    sll_real64,dimension(:),pointer :: buf,Xstar,node_pos,coeffs
    sll_int32,dimension(:),pointer :: ibuf 
    
    
    dx = 1._f64/real(N,f64)
    
    allocate(buf(10*N))
    allocate(ibuf(N))
    allocate(node_pos(-2:N+2),coeffs(-2:N+2))
    allocate(Xstar(1:N+1))
    
    
    
    node_pos(0:N)=node_positions(1:N+1)
    
    
    !do i=1,N+1
    !  print *,i,node_positions(i)
    !enddo

    !do i=1,N+1
    !  print *,i,Xstar(i),f(i)
    !enddo
    
    !print *,dx
    
    do i=1,N+1    
      Xstar(i) = node_positions(i)-alpha
    enddo
    
    
    do i=1,N+1
      do while (Xstar(i).gt.1._f64)
        Xstar(i) = Xstar(i)-1._f64
      end do
      do while (Xstar(i).lt.0._f64)
        Xstar(i) = Xstar(i)+1._f64
      end do    
    enddo

   

    !from f compute the mean
!    do i=1,N
!      f(i)=f(i)*(node_positions(i+1)-node_positions(i))/dx
!    enddo
    
    
    !we compute the splines coefficients by solving the LU decomposition
!    M=0._f64
!    do i=1,N
!      M=M+f(i)
!    enddo
!    !M=M/real(N,rk)
!    do i=1,N
!      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
!    enddo    
!    tmp=f(1)
!    f(1)=0._f64
!    do i=2,N
!      tmp2=f(i)
!      f(i)=f(i-1)+tmp
!      tmp=tmp2
!    enddo
    
    call setup_spline_nonunif_1D_periodic_aux( node_pos, N, buf, ibuf)
    call compute_spline_nonunif_1D_periodic_aux( f, N, buf, ibuf, coeffs )
    call interpolate_array_value_nonunif_aux( Xstar, f, N, node_pos, coeffs,N)
    
    
!    tmp=f(1)
!    do i=1,N-1
!      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
!    enddo
!    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
!    do i=1,N
!      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
!    enddo

    f(N+1) = f(1)
    
    deallocate(buf)
    deallocate(ibuf)
    deallocate(node_pos,coeffs)
    deallocate(Xstar)
    
  end subroutine constant_advection_spl_non_unif_per
  
  subroutine function_to_primitive_old(f,node_positions,N,M)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    sll_int32,intent(in):: N
    sll_real64,intent(out)::M
    sll_int32::i
    sll_real64::dx,tmp,tmp2
    dx = 1._f64/real(N,f64)
    
    !from f compute the mean
    do i=1,N
      f(i)=f(i)*(node_positions(i+1)-node_positions(i))/dx
    enddo
    
    
    M=0._f64
    do i=1,N
      M=M+f(i)
    enddo
    !M=M/real(N,rk)
    do i=1,N
      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
    enddo    
    tmp=f(1)
    f(1)=0._f64
    do i=2,N
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo



  end subroutine function_to_primitive_old




  subroutine primitive_to_function_old(f,node_positions,N,M)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    sll_int32,intent(in):: N
    sll_real64,intent(in)::M
    sll_int32::i
    sll_real64::dx,tmp,tmp2
    dx = 1._f64/real(N,f64)
    
    tmp=f(1)
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
    do i=1,N
      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
    enddo

    f(N+1) = f(1)



  end subroutine primitive_to_function_old



  subroutine function_to_primitive(f,node_positions,N,M)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    sll_int32,intent(in):: N
    sll_real64,intent(out)::M
    sll_int32::i
    sll_real64::dx,tmp,tmp2
    dx = 1._f64/real(N,f64)
        
    !from f compute the mean
    M=0._f64
    do i=1,N
      M=M+f(i)*(node_positions(i+1)-node_positions(i))
    enddo
    
    f(1)=(f(1)-M)*(node_positions(2)-node_positions(1))
    tmp=f(1)
    f(1)=0._f64
    do i=2,N+1
      f(i)=(f(i)-M)*(node_positions(i+1)-node_positions(i))
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo    
    
    
    
    !print *,f(1),f(N+1) 

  end subroutine function_to_primitive




  subroutine primitive_to_function(f,node_positions,N,M)
    implicit none
    sll_real64,dimension(:),allocatable::f,node_positions
    sll_int32,intent(in):: N
    sll_real64,intent(in)::M
    sll_int32::i
    sll_real64::dx,tmp,tmp2
    dx = 1._f64/real(N,f64)
    
    tmp=f(1)
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
    do i=1,N
      f(i)=f(i)/(node_positions(i+1)-node_positions(i))
    enddo

    f(N+1) = f(1)



  end subroutine primitive_to_function








  subroutine csl_constant_advection_spl_non_unif_per(f,alpha,node_positions,N)
    !alpha and node_positions are normalized to [0,1]
    !use numeric_constants
    !use cubic_non_uniform_splines
    implicit none
    
    sll_real64,dimension(:),allocatable::f,node_positions
    !type(cubic_nonunif_spline_1D), pointer :: spl_per
    sll_int32,intent(in):: N
    sll_real64,intent(in)::alpha
    sll_real64 :: dx
    sll_int32  :: i
    sll_real64 :: M,tmp,tmp2
    !temporary allocations
    sll_real64,dimension(:),pointer :: buf,Xstar,node_pos,coeffs
    sll_int32,dimension(:),pointer :: ibuf 
    
    
    dx = 1._f64/real(N,f64)
    
    allocate(buf(10*N))
    allocate(ibuf(N))
    allocate(node_pos(-2:N+2),coeffs(-2:N+2))
    allocate(Xstar(1:N+1))
    

    

    
    
    node_pos(0:N)=node_positions(1:N+1)
    
    
    !do i=1,N+1
    !  print *,i,node_positions(i)
    !enddo

    !do i=1,N+1
    !  print *,i,Xstar(i),f(i)
    !enddo
    
    !print *,dx
    
    do i=1,N+1    
      Xstar(i) = node_positions(i)-alpha
    enddo
    
    
    do i=1,N+1
      do while (Xstar(i).gt.1._f64)
        Xstar(i) = Xstar(i)-1._f64
      end do
      do while (Xstar(i).lt.0._f64)
        Xstar(i) = Xstar(i)+1._f64
      end do    
    enddo

   

    !from f compute the mean
    do i=1,N
      f(i)=f(i)*(node_positions(i+1)-node_positions(i))/dx
    enddo
    
    
    !we compute the splines coefficients by solving the LU decomposition
    M=0._f64
    do i=1,N
      M=M+f(i)
    enddo
    !M=M/real(N,rk)
    do i=1,N
      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
    enddo    
    tmp=f(1)
    f(1)=0._f64
    do i=2,N
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo

    !do i=1,N+1
    !  print *,i,Xstar(i),node_pos(i-1),Xstar(i)-node_pos(i-1)
    !enddo
    !stop


    
    call setup_spline_nonunif_1D_periodic_aux( node_pos, N, buf, ibuf)
    call compute_spline_nonunif_1D_periodic_aux( f, N, buf, ibuf, coeffs )
    call interpolate_array_value_nonunif_aux( Xstar(1:N), f(1:N), N, node_pos, coeffs,N)
    
    
    tmp=f(1)
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
    do i=1,N
      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
    enddo

    f(N+1) = f(1)
    
    
    deallocate(buf)
    deallocate(ibuf)
    deallocate(node_pos,coeffs)
    deallocate(Xstar)
    
  end subroutine csl_constant_advection_spl_non_unif_per


  subroutine interpolate_array_non_unif_per(f_new,new_positions,N_new,f,node_positions,N)
    !Xstar and node_positions are normalized to [0,1]
    !use numeric_constants
    !use cubic_non_uniform_splines
    implicit none
    
    sll_real64,dimension(:),allocatable::f,node_positions,new_positions,f_new
    !type(cubic_nonunif_spline_1D), pointer :: spl_per
    sll_int32,intent(in):: N,N_new
    sll_real64 :: dx
    sll_int32  :: i
    sll_real64 :: M,tmp,tmp2
    !temporary allocations
    sll_real64,dimension(:),pointer :: buf,Xstar,node_pos,coeffs
    sll_int32,dimension(:),pointer :: ibuf 
    
    
    dx = 1._f64/real(N,f64)
    
    allocate(buf(10*N))
    allocate(ibuf(N))
    allocate(node_pos(-1:N+2),coeffs(-1:N+2))
    allocate(Xstar(1:N+1))
    
    
    
    node_pos(0:N)=node_positions(1:N+1)
    
    
    !do i=1,N+1
    !  print *,i,node_positions(i)
    !enddo

    !do i=1,N+1
    !  print *,i,Xstar(i),f(i)
    !enddo
    
    !print *,dx
    
    !do i=1,N+1    
    !  Xstar(i) = node_positions(i)-alpha
    !enddo
    
    
    !do i=1,N+1
    !  do while (Xstar(i).gt.1._f64)
    !    Xstar(i) = Xstar(i)-1._f64
    !  end do
    !  do while (Xstar(i).lt.0._f64)
    !    Xstar(i) = Xstar(i)+1._f64
    !  end do    
    !enddo

   

    !from f compute the mean
!    do i=1,N
!      f(i)=f(i)*(node_positions(i+1)-node_positions(i))/dx
!    enddo
    
    
    !we compute the splines coefficients by solving the LU decomposition
!    M=0._f64
!    do i=1,N
!      M=M+f(i)
!    enddo
!    !M=M/real(N,rk)
!    do i=1,N
!      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
!    enddo    
!    tmp=f(1)
!    f(1)=0._f64
!    do i=2,N
!      tmp2=f(i)
!      f(i)=f(i-1)+tmp
!      tmp=tmp2
!    enddo
    
    call setup_spline_nonunif_1D_periodic_aux( node_pos, N, buf, ibuf)
    call compute_spline_nonunif_1D_periodic_aux( f, N, buf, ibuf, coeffs )
    call interpolate_array_value_nonunif_aux( new_positions, f_new, N_new, node_pos, coeffs,N)
    
    
!    tmp=f(1)
!    do i=1,N-1
!      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
!    enddo
!    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
!    do i=1,N
!      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
!    enddo

    !f(N+1) = f(1)
    
    deallocate(buf)
    deallocate(ibuf)
    deallocate(node_pos,coeffs)
    deallocate(Xstar)
    
  end subroutine interpolate_array_non_unif_per





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! functions for computing integral
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  function compute_non_unif_integral(integration_points,N_points,rho_case)
    real(f64) :: compute_non_unif_integral
    real(f64),dimension(:,:),allocatable :: integration_points
    integer,intent(in) :: N_points,rho_case
    if(rho_case==0)then
      compute_non_unif_integral= compute_non_unif_integral_rectangle(integration_points,N_points)
    endif
    if(rho_case==1)then
      compute_non_unif_integral= compute_non_unif_integral_trapezoid(integration_points,N_points)
    endif
    if(rho_case==2)then
      !for conservative on dual mesh (x_{i+1/2}=(x_i+x_{i+1})/2)
      compute_non_unif_integral=compute_non_unif_integral_trapezoid(integration_points,N_points)
      compute_non_unif_integral=compute_non_unif_integral&
      &-0.5_f64*integration_points(2,N_points)*(integration_points(1,2)-integration_points(1,1))&
      &+0.5_f64*integration_points(2,1)*(integration_points(1,N_points)-integration_points(1,N_points-1))
    endif
!    if(rho_case==3)then
!      compute_non_unif_integral=compute_non_unif_integral_gaussian(integration_points,N_points)
!    endif      
!    if(rho_case==4)then
!      compute_non_unif_integral=compute_non_unif_integral_gaussian_sym(integration_points,N_points)
!    endif        
!    if(rho_case==5)then
!      compute_non_unif_integral=compute_non_unif_integral_spline_per(integration_points,N_points)
!    endif
  
end  function compute_non_unif_integral

!function compute_non_unif_integral_hermite(integration_points,N_points,d)
!  real(f64) :: compute_non_unif_integral_trapezoid
!  real(f64),dimension(:,:),allocatable :: integration_points
!  integer,intent(in) :: N_points,d
!  
!end function compute_non_unif_integral_hermite
!

function compute_non_unif_integral_trapezoid(integration_points,N_points)
  real(f64) :: compute_non_unif_integral_trapezoid
  real(f64),dimension(:,:),allocatable :: integration_points
  integer,intent(in) :: N_points
  integer :: i
  real(f64) :: tmp,x1,x2,fval1,fval2
  compute_non_unif_integral_trapezoid = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    if(x2<x1)then
      print *,i,'bad integration points x1=',x1,'x2=',x2
      stop
    endif
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    tmp = 0.5_f64*(fval1+fval2)*(x2-x1)
    compute_non_unif_integral_trapezoid=compute_non_unif_integral_trapezoid+tmp
  enddo
  
  
end  function compute_non_unif_integral_trapezoid


function compute_non_unif_integral_rectangle(integration_points,N_points)
  real(f64) :: compute_non_unif_integral_rectangle
  real(f64),dimension(:,:),allocatable :: integration_points
  integer,intent(in) :: N_points
  integer :: i
  real(f64) :: tmp,x1,x2,fval
  compute_non_unif_integral_rectangle = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    if(x2<x1)then
      print *,i,'bad integration points x1=',x1,'x2=',x2
      stop
    endif
    fval = integration_points(2,i)
    !fval2 = integration_points(2,i+1)
    tmp = fval*(x2-x1)
    compute_non_unif_integral_rectangle=compute_non_unif_integral_rectangle+tmp
  enddo
  
  
end  function compute_non_unif_integral_rectangle



  subroutine poisson1dper_init(coefd,N)  
    integer,intent(in)::N
    real(f64),dimension(2*N+15),intent(out)::coefd
    call  dffti(N,coefd)
  end subroutine poisson1dper_init




  subroutine computerhoper(ftab,dom,Nx,Nv,Nz,rho)
    integer,intent(in)::Nx,Nv,Nz
    real(f64),dimension(0:Nx-1,0:Nv-1),intent(in)::ftab
    real(f64),dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:Nz-1),intent(out)::rho
    real(f64)::dv
    integer::i,j
    dv=dom(1,1)/real(Nv,f64)
    rho=0._f64
    do j=0,Nv-1
      do i=0,Nx-1
        rho(i)=rho(i)+ftab(i,j)
      enddo
    enddo  
    rho=rho*dv
  end subroutine computerhoper



  subroutine poisson1dper(coefd,E,L,N)
    integer,intent(in)::N
    real(f64),dimension(2*N+15),intent(in)::coefd
    real(f64),dimension(N),intent(inout)::E
    real(f64),intent(in)::L
    integer::i
    real(f64)::re,im,tmp
    call dfftf(N,E,coefd)
    tmp=0.5_f64*(L/sll_pi)/real(N,f64);
        
    E(1)=0._f64
    do i=1,(N-2)/2
      re=E(2*i);im=E(2*i+1)
      E(2*i)=tmp/real(i,f64)*im
      E(2*i+1)=-tmp/real(i,f64)*re
    enddo
    if(mod(N,2)==0)E(N)=0._f64
    call dfftb(N,E,coefd)
  
  end subroutine poisson1dper


  subroutine printgp1dper(xmin,xmax,ftab,N)
    integer,intent(in)::N
    real(f64),intent(in)::xmin,xmax
    real(f64),dimension(0:N-1),intent(in)::ftab
    integer::i
    real(f64)::x,dx
    open(unit=900,file='f.dat')
    dx=(xmax-xmin)/real(N,f64)    
    do i=0,N-1
      x=xmin+real(i,f64)*dx      
      write(900,*) x,ftab(i)
    enddo
    write(900,*) xmax,ftab(0)
    close(900)  
  end subroutine printgp1dper

  subroutine printgp1dnuper(X,ftab,N)
    integer,intent(in)::N
    real(f64),dimension(0:N)::X
    real(f64),dimension(0:N-1),intent(in)::ftab
    integer::i
    open(unit=900,file='f.dat')
       
    do i=0,N-1
      write(900,*) X(i),ftab(i)
    enddo
    write(900,*) X(N),ftab(0)
    close(900)  
  end subroutine printgp1dnuper

  subroutine printgp1dper_non_unif(X,ftab,N,step,filename)
    integer,intent(in)::N,step
    real(f64),dimension(0:N)::X
    real(f64),dimension(0:N-1),intent(in)::ftab
    integer::i
    character(len=*),intent(in)::filename
    character*80,str,str2
    
    write(str2,*) step
    !write(str,*) trim(adjustl(filename//trim(adjustl((str2)))//'.dat'))
    open(unit=900,file='f.dat')       
    do i=0,N-1
      write(900,*) X(i),ftab(i)
    enddo
    write(900,*) X(N),ftab(0)
    close(900)
    write(str,*) 'mv f.dat '//filename//trim(adjustl((str2)))//'.dat';call system(str)  
  end subroutine printgp1dper_non_unif



  subroutine print2dper(ftab,Nx,Ny,visucase,step,filename)
    integer,intent(in)::Nx,Ny,visucase,step
    !real(f64),dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    character(len=*),intent(in)::filename
    character*80,str,str2
    
        
    if(visucase==0)then
      !gnuplot
      call printgp2dper(ftab,Nx,Ny)
      write(str2,*) step
      !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
      write(str,*) 'mv f.dat '//filename//trim(adjustl((str2)))//'.dat';call system(str)
    endif
    if(visucase==1)then
      !vtk
      call printvtk2dper(ftab,Nx,Ny)
      write(str2,*) step
      write(str,*) 'mv f.vtk '//filename//trim(adjustl((str2)))//'.vtk';call system(str)
    endif
  end subroutine print2dper

  subroutine printgp2dper(ftab,Nx,Ny)
    integer,intent(in)::Nx,Ny
    !real(f64),dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    integer::i,j
    real(f64)::z(0:1),dz(0:1)
    dz(0)=1._f64/real(Nx,f64);dz(1)=1._f64/real(Ny,f64)
    open(unit=900,file='f.dat')
    do j=0,Ny-1
      do i=0,Nx-1
        z(0)=real(i,f64)*dz(0)
        z(1)=real(j,f64)*dz(1)
        write(900,*) z(0),z(1),ftab(i,j)
      enddo
      i=Nx
      z(0)=real(i,f64)*dz(0)
      z(1)=real(j,f64)*dz(1)
      write(900,*) z(0),z(1),ftab(0,j)      
      write(900,*) ''      
    enddo
    j=Ny
    do i=0,Nx-1
      z(0)=real(i,f64)*dz(0)
      z(1)=real(j,f64)*dz(1)
      write(900,*) z(0),z(1),ftab(i,0)
    enddo
    i=Nx
    z(0)=real(i,f64)*dz(0)
    z(1)=real(j,f64)*dz(1)
    write(900,*) z(0),z(1),ftab(0,0)	  
    write(900,*) ''	       
    close(900)  
  end subroutine printgp2dper

  subroutine printvtk2dper(ftab,Nx,Ny)
    integer,intent(in)::Nx,Ny
    !real(f64),dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    integer::i,j
    real(f64)::z(0:1),dz(0:1)
    dz(0)=1._f64/real(Nx,f64);dz(1)=1._f64/real(Ny,f64)
    !open(unit=900,file='f.vtk')
    open(unit=900,file='f.vtk',form='formatted')
    write(900,'(A)')                  '# vtk DataFile Version 2.0'
    write(900,'(A)')                  'Exemple'
    write(900,'(A)')                  'ASCII'
    write(900,'(A)')                  'DATASET STRUCTURED_POINTS'
    write(900,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx+1,' ', Ny+1,' ', 1
    write(900,'(A,I0,A,I0,A,I0)') 'ORIGIN ', 0,' ' , 0,' ' , 0
    !write(900,'(A,F10.4,A,F10.4,A,F10.4)') 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*) 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*)
    write(900,'(A,I0)')           'POINT_DATA ',(Nx+1)*(Ny+1)
    write(900,'(A,I0)')           'SCALARS f float ',1
    write(900,'(A)')                  'LOOKUP_TABLE default'
    
    do j=0,Ny-1
      do i=0,Nx-1
        z(0)=real(i,f64)*dz(0)
        z(1)=real(j,f64)*dz(1)
        !write(900,'(F0.8)') ftab(i,j)
        write(900,*) ftab(i,j)
      enddo
      i=Nx
      z(0)=real(i,f64)*dz(0)
      z(1)=real(j,f64)*dz(1)
      !write(900,'(F0.8)') ftab(0,j)            
      write(900,*) ftab(0,j)            
    enddo
    j=Ny
    do i=0,Nx-1
      z(0)=real(i,f64)*dz(0)
      z(1)=real(j,f64)*dz(1)
      !write(900,'(F0.8)') ftab(i,0)
      write(900,*) ftab(i,0)
    enddo
    i=Nx
    z(0)=real(i,f64)*dz(0)
    z(1)=real(j,f64)*dz(1)
    !write(900,'(F0.8)') ftab(0,0)	  	       
    write(900,*) ftab(0,0)	  	       
    close(900)  
  end subroutine printvtk2dper

  subroutine constant_advection_size(sizebuf,interp,N)
    integer,intent(in)::N
    integer,intent(out)::sizebuf
    integer,dimension(2),intent(in)::interp
    sizebuf=1
         
    if(interp(1)==2)then
      sizebuf=3*N+15
    endif
    if(interp(1)==1)then
      sizebuf=N
      if(interp(2)==3)then
        sizebuf=5*N
      endif
    endif
    if(interp(1)==0)then
      sizebuf=N
    endif
  end subroutine constant_advection_size
  

  subroutine constant_advection_init(buf,sizebuf,N,interp)
    integer,intent(in)::sizebuf,N
    integer,dimension(2)::interp    
    real(f64),dimension(0:sizebuf-1),intent(out)::buf
    
    if(interp(1)==2)then
      call fourier1dperlagodd_init(buf,sizebuf,N)
    endif
    if(interp(1)==1)then
      call hermiteper1d_init(buf,sizebuf,interp(2),N)
    endif
    if(interp(1)==0)then
      call hermite_finite_difference_per_init(buf,sizebuf,N,interp(2))
    endif
    
  end subroutine constant_advection_init

  subroutine constant_advection_solve(buf,sizebuf,f,N,alpha,interp)
    integer,intent(in)::sizebuf,N
    integer,dimension(2)::interp    
    real(f64),dimension(0:sizebuf-1),intent(inout)::buf
    real(f64),dimension(1:N),intent(inout)::f
    real(f64),intent(in)::alpha
    
    if(interp(1)==2)then
      call fourier1dperlagodd(buf,sizebuf,f,N,alpha,(interp(2)-1)/2)
    endif
    if(interp(1)==1)then
      call hermiteper1d(buf,sizebuf,f,N,alpha,interp(2))
    endif
    if(interp(1)==0)then
      call hermite_finite_difference_per(buf,sizebuf,f,N,alpha,interp(2))
    endif        

  end subroutine constant_advection_solve

  subroutine hermite_finite_difference_per_init(buf,sizebuf,N,p)
    integer,intent(in)::sizebuf,N,p    
    real(f64),dimension(0:sizebuf-1),intent(inout)::buf
    
  end subroutine hermite_finite_difference_per_init

  subroutine compute_finite_difference_per(f,df,N,w,p)
    integer,intent(in)::N,p    
    real(f64),dimension(-p/2:(p+1)/2),intent(in)::w
    real(f64),dimension(0:N-1),intent(in)::f
    real(f64),dimension(2,0:N-1),intent(out)::df
    integer::r,s,i,ii
    real(f64)::tmp
    r=-p/2
    s=(p+1)/2


    do i=0,N-1
      tmp=0._f64
      do ii=r,s
        tmp=tmp+w(ii)*f(modulo(i+ii+N,N))
      enddo
      df(1,i)=tmp
      if(-r/=s)then
        tmp=0._f64
        do ii=r,s
          tmp=tmp-w(r+s-ii)*f(modulo(i+ii+N,N))
        enddo
        df(2,i)=tmp
      endif        
    enddo
    if(-r==s)then      
      do i=0,N-2
        df(2,i)=df(1,i+1)
      enddo
      df(2,N-1)=df(1,0)
    endif



    
  end subroutine compute_finite_difference_per

  subroutine compute_finite_difference_init(w,p)
    integer,intent(in)::p    
    real(f64),dimension(-p/2:(p+1)/2),intent(out)::w
    integer::r,s,i,j
    real(f64)::tmp

    r=-p/2
    s=(p+1)/2

!    if(modulo(p,2)==0)then
!      r=-p/2
!      s=p/2
!    else
!      r=-(p-1)/2
!      s=(p+1)/2
!    endif
        
    
    !maple code for generation of w
    !for k from r to -1 do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
    !od:
    !for k from 1 to s do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
    !od:
    !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):
    
    do i=r,-1
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

    do i=1,s
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo
    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp
    
  end subroutine compute_finite_difference_init

  subroutine hermite_finite_difference_per(buf,sizebuf,f,N,alpha,p)
    integer,intent(in)::sizebuf,N,p    
    real(f64),dimension(0:sizebuf-1),intent(inout)::buf
    real(f64),dimension(0:N-1),intent(inout)::f
    real(f64),intent(in)::alpha
    real(f64),dimension(:),allocatable::w
    real(f64),dimension(:,:),allocatable::df
    integer::i,r,s,ii,i0,i1
    real(f64)::tmp,x,ww(0:1,0:1)

    r=-p/2
    s=(p+1)/2
    allocate(w(r:s))
    allocate(df(2,0:N-1))
    call compute_finite_difference_init(w,p)
!    print *,'w',w
    !do ii=r,s
    !  print *,ii,r+s-ii,w(ii),-w(r+s-ii)
    !enddo
    
    !stop
    
    !we compute the derivatives
    buf(0:N-1)=f(0:N-1)
    call compute_finite_difference_per(buf(0:N-1),df(1:2,0:N-1),N,w,p)
    !we make the interpolation on uniform grid with derivative info
    !localization
    x=-alpha
    !x=0._f64
    x=x-real(floor(x),f64)
    x=x*real(N,f64)
    i0=floor(x)
    if(i0==N)then
      x=0._f64;i0=0
    endif
    x=x-real(i0,f64)
    ww(1,0)=x**2*(3._f64-2._f64*x)
    ww(0,0)=1._f64-ww(1,0)
    ww(0,1)=x*(1._f64-x)**2
    ww(1,1)=x**2*(x-1._f64)
    !update f from buf and df
    i1=i0
    do i=0,N-1 
      i1=i0+1
      if(i1>=N)then
        i1=i1-N
      endif        
      !we use buf(i0),buf(i1),df(1,i0),df(2,i0) 
      f(i)=ww(0,0)*buf(i0)+ww(1,0)*buf(i1)+ww(0,1)*df(1,i0)+ww(1,1)*df(2,i0)
      i0=i1
    enddo
    deallocate(w)
    
  end subroutine hermite_finite_difference_per

  subroutine fourier1dperlagodd_init(buf,sizebuf,N)  
    integer,intent(in)::N,sizebuf
    real(f64),dimension(0:sizebuf-1)::buf
    call dffti(N,buf(N:3*N+14))
  end subroutine fourier1dperlagodd_init


  subroutine fourier1dperlagodd(buf,sizebuf,E,N,alpha,d)
    integer,intent(in)::N,sizebuf,d
    real(f64),dimension(0:sizebuf-1),intent(inout)::buf
    real(f64),dimension(1:N),intent(inout)::E
    real(f64),intent(in)::alpha
    integer::i,ix
    real(f64)::rea,ima,reb,imb,tmp,x,a
    !localization
    x=alpha
    !if(abs(x)<1.e-15_f64)then
    !  print *,x,N
    !  x=x-real(floor(x),f64)
    !  print *,x
    !  x=x*real(N,f64)
    !  print *,x
    !  ix=floor(x)
    !  x=x-real(ix,f64)    
    !  print *,x,ix
    !  x=0._f64
    !endif
    x=x-real(floor(x),f64)
    !if(x==1)then
    !  x=0._f64
    !endif  
    x=x*real(N,f64)
    ix=floor(x)
    if(ix==N)then
      x=0._f64;ix=0
    endif
    x=x-real(ix,f64)    
    do i=0,N-1
      buf(i)=0._f64
    enddo

    a=1._f64;
    do i=2,d
      a=a*(x*x-real(i,f64)*real(i,f64))/(real(d,f64)*real(d,f64))
    enddo
    a=a*(x+1._f64)/real(d,f64)
    a=a*(x-real(d,f64)-1._f64)/real(d,f64)
    buf(ix)=a*(x-1._f64)/real(d,f64)
    buf(mod(ix+1,N))=a*x/real(d,f64)
    a=a*x*(x-1._f64)/(real(d,f64)*real(d,f64))  
    do i=-d,-1
      buf(mod(i+ix+N,N))=a/((x-real(i,f64))/real(d,f64))
    enddo  
    do i=2,d+1
      buf(mod(i+ix+N,N))=a/((x-real(i,f64))/real(d,f64));
    enddo
    a=1._f64;
    do i=-d,d+1
      buf(mod(i+ix+N,N))=buf(mod(i+ix+N,N))*a
      a=a*real(d,f64)/real(d+i+1,f64)
    enddo
    a=1._f64;
    do i=d+1,-d,-1
      buf(mod(i+ix+N,N))=buf(mod(i+ix+N,N))*a
      a=a*real(d,f64)/real(i-1-d-1,f64)
    enddo


    call dfftf(N,buf(0:N-1),buf(N:3*N+14))

    call dfftf(N,E,buf(N:3*N+14))
    tmp=1._f64/real(N,f64);            
    E(1)=E(1)*tmp*buf(0)
    do i=1,(N-2)/2
      rea=E(2*i);ima=E(2*i+1)
      reb=tmp*buf(2*i-1);imb=tmp*buf(2*i);
      E(2*i)=rea*reb-ima*imb
      E(2*i+1)=rea*imb+reb*ima
    enddo
    if(mod(N,2)==0)E(N)=E(N)*tmp*buf(N-1)
    call dfftb(N,E,buf(N:3*N+14))
  end subroutine fourier1dperlagodd

  subroutine hermiteper1d_init(buf,sizebuf,interp,N)
    integer,intent(in)::sizebuf,interp,N    
    real(f64),dimension(0:sizebuf-1),intent(out)::buf
    if(interp==3)then
      !call splcoefper1d0(buf(2*N:5*N-1),N)
      !call splcoefper1d0new(buf(0:3*N-1),N)
      call splcoefper1d0new(buf(2*N:5*N-1),N)
    endif
  end subroutine hermiteper1d_init



  subroutine hermiteper1d(buf,sizebuf,p,N,alpha,interp)
    integer,intent(in)::N,interp,sizebuf
    real(f64),intent(in)::alpha
    real(f64),dimension(0:N-1),intent(inout)::p
    real(f64),dimension(0:sizebuf-1),intent(inout)::buf
    real(f64)::x,w(0:2),fbar,df0,df1,tmp
    integer::i,ix,i0,i1,i2,im1,im2,im3
    
    !localization
    x=-alpha
    !x=0._f64
    x=x-real(floor(x),f64)
    x=x*real(N,f64)
    ix=floor(x)
    if(ix==N)then
      x=0._f64;ix=0
    endif
    x=x-real(ix,f64)
    !print *,x,ix
    
    w(0)=x*(1._f64-x)*(1._f64-x)
    w(1)=x*x*(x-1._f64)
    w(2)=x*x*(3._f64-2._f64*x)
    i0=ix
    i1=mod(ix+1,N)
    i2=mod(ix+2,N)
    im1=mod(ix+N-1,N)
    im2=mod(ix+N-2,N)
    im3=mod(ix+N-3,N)
    
    !do i=0,N-1,
    !  buf(i)=p(i)
    !enddo
    buf(0:N-1)=p(0:N-1)
    
    select case (interp)
      case(0)!pfc without limiter
      !w0:=x*(1-x)**2:w1:=x**2*(x-1):w2:=x**2*(3-2*x):
      !df0:=5/6*f1-f2/6+f0/3:df1:=5/6*f1+f2/3-f0/6:
      !collect(sort(expand(w0*df0+w1*df1+w2*f1),x),x);
      !B0:=x*f1-x*(x-2)*(x-1)/6*(f1-f0)+x*(x-1)*(x+1)/6*(f2-f1);
      !B:=collect(sort(expand(B0,x)),x);
        do i=0,N-1
	  fbar=buf(i0)
	  df0=(5._f64/6._f64)*buf(i0)-(1._f64/6._f64)*buf(i1)+(2._f64/6._f64)*buf(im1)
	  df1=(5._f64/6._f64)*buf(i0)+(2._f64/6._f64)*buf(i1)-(1._f64/6._f64)*buf(im1)
	  p(i)=w(0)*df0+w(1)*df1+w(2)*fbar;
	  im1=i0;i0=i1;i1=i1+1;if(i1>=N)i1=i1-N;
	enddo
      case(1)!ppm without limiter
        tmp=(7._f64/12._f64)*(buf(i0)+buf(im1))-(1._f64/12._f64)*(buf(i1)+buf(im2))
	do i=0,N-1
	  fbar=buf(i0)
	  df0=tmp
	  im2=im1;im1=i0;i0=i1;i1=i1+1;if(i1>=N)i1=i1-N;
	  df1=(7._f64/12._f64)*(buf(i0)+buf(im1))-(1._f64/12._f64)*(buf(i1)+buf(im2))
	  tmp=df1
	  p(i)=w(0)*df0+w(1)*df1+w(2)*fbar;
	enddo
      case(2)!ppm2 without limiter
        tmp=(37._f64/60._f64)*(buf(i0)+buf(im1))-(8._f64/60._f64)*(buf(i1)+buf(im2))+(1._f64/60._f64)*(buf(i2)+buf(im3))
	do i=0,N-1
	  fbar=buf(i0)
	  df0=tmp
	  im3=im2;im2=im1;im1=i0;i0=i1;i1=i2;i2=i2+1;if(i2>=N)i2=i2-N;
	  df1=(37._f64/60._f64)*(buf(i0)+buf(im1))-(8._f64/60._f64)*(buf(i1)+buf(im2))+(1._f64/60._f64)*(buf(i2)+buf(im3))
	  tmp=df1
	  p(i)=w(0)*df0+w(1)*df1+w(2)*fbar;
	enddo
      case(3)!psm without limiter
        !buf(N:2*N-1)=df(0:N-1) computed as splines coefficients of 0.5*(f_i+f_{i-1}), i=0..N-1
	buf(N+0)=0.5_f64*(buf(0)+buf(N-1))
        do i=1,N-1
          buf(N+i)=0.5_f64*(buf(i)+buf(i-1))
        enddo
        !call splcoefper1d(buf(N:2*N-1),buf(2*N:5*N-1),N)
        call splcoefper1dnew(buf(N:2*N-1),buf(2*N:5*N-1),N)
	!print *,buf(N:2*N-1)
	tmp=buf(N+i0)
	do i=0,N-1
	  fbar=buf(i0)
	  df0=tmp
	  i0=i0+1;if(i0>=N)i0=i0-N;
	  df1=buf(N+i0)
	  tmp=df1
	  p(i)=w(0)*df0+w(1)*df1+w(2)*fbar;
	enddo
      case default
    end select

    i0=ix;
    fbar=p(0);
    do i=0,N-2
      p(i)=buf(i0)+p(i+1)-p(i)
      i0=i0+1;if(i0>=N)i0=i0-N
    enddo
    p(N-1)=buf(i0)+fbar-p(N-1)
        
  end subroutine hermiteper1d


  subroutine splcoefper1d0new(luper,N)
    integer,intent(in)::N
    real(f64),dimension(0:3*N-1),intent(out)::luper
    integer::i
    
    luper(0+3*0)=4._f64
    luper(2+3*0)=0.25_f64
    do i=0,N-2
      luper(1+3*i)=1._f64/luper(0+3*i)
      luper(0+3*(i+1))=4._f64-luper(1+3*i)
      luper(2+3*(i+1))=-luper(2+3*i)/luper(0+3*(i+1))
    enddo
    luper(0+3*(N-1))=luper(0+3*(N-1))-(luper(1+3*(N-2))+2._f64*luper(2+3*(N-2)))  
    do i=0,N-1
      luper(0+3*i)=1._f64/luper(0+3*i)
    enddo
  end subroutine splcoefper1d0new



  subroutine splcoefper1dnew(f,luper,N)
    integer,intent(in)::N
    real(f64),dimension(0:3*N-1),intent(in)::luper
    real(f64),dimension(0:N-1),intent(inout)::f
    integer::i
    do i=0,N-1;f(i)=6._f64*f(i);enddo;
    do i=1,N-1
      f(i)=f(i)-f(i-1)*luper(1+3*(i-1))
    enddo
    do i=0,N-2
      f(N-1)=f(N-1)-luper(2+3*i)*f(i)
    enddo
    f(N-1)=f(N-1)*luper(0+3*(N-1));f(N-2)=luper(0+3*(N-2))*(f(N-2)-(1._f64-luper(2+3*(N-3)))*f(N-1))
    do i=N-3,1,-1
      f(i)=luper(0+3*i)*(f(i)-f(i+1)+luper(2+3*(i-1))*f(N-1))
    enddo
    f(0)=luper(0+3*0)*(f(0)-f(1)-f(N-1));
  end subroutine splcoefper1dnew



  subroutine thdiagper(ftab,E,dom,Nx,Nv,Nz,nbdiag,thf)
    integer,intent(in)::Nx,Nv,Nz,nbdiag
    real(f64),dimension(0:Nx-1,0:Nv-1),intent(in)::ftab
    real(f64),dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:Nz-1),intent(in)::E
    real(f64),dimension(0:nbdiag-1),intent(out)::thf
    real(f64)::dx,dv,tmp,l1,l2,ee,mass,ke,mv,v,x,et,rms
    integer::i,j
    dx=dom(1,0)/real(Nx,f64);dv=dom(1,1)/real(Nv,f64)
    
    l1=0._f64;l2=0._f64;ee=0._f64;mass=0._f64;ke=0._f64;mv=0._f64;rms=0._f64
    do j=0,Nv-1
      do i=0,Nx-1
    	v=dom(0,1)+real(j,f64)*dv
	x=dom(0,0)+real(i,f64)*dx
	l1=l1+abs(ftab(i,j))
    	l2=l2+abs(ftab(i,j))*abs(ftab(i,j))
    	mass=mass+ftab(i,j)
	ke=ke+v*v*ftab(i,j)
	mv=mv+v*ftab(i,j)      
        rms=rms+x*x*ftab(i,j)
      enddo
    enddo  
    
    do i=0,Nx-1
      ee=ee+abs(E(i))*abs(E(i))
    enddo	
    l1=l1*dx*dv;l2=sqrt(l2*dx*dv);mass=mass*dx*dv;ee=sqrt(ee*dx)
    ke=ke*dx*dv
    et=ke+ee*ee;mv=mv*dx*dv;rms=sqrt(rms*dx*dv)

   thf(0)=ee;thf(1)=mass;thf(2)=l1;thf(3)=l2;thf(4)=et;thf(5)=mv;thf(6)=ke;thf(7)=rms
    
  end subroutine thdiagper



  subroutine PFenvelope(S, t, tflat, tL, tR, twL, twR, t0, &
       turn_drive_off)

    ! DESCRIPTION
    ! -----------
    ! S: the wave form at a given point in time. This wave form is 
    !    not scaled (its maximum value is 1).
    ! t: the time at which the envelope is being evaluated
    ! tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
    !    envelope, defined in the main portion of this program.
    ! turn_drive_off: 1 if the drive should be turned off after a time
    !    tflat, and 0 otherwise

    real(f64), intent(in) :: t, tflat, tL, tR, twL, twR, t0
    real(f64), intent(out) :: S
    logical, intent(in) :: turn_drive_off
    ! local variables
    integer :: i 
    real(f64) :: epsilon

    ! The envelope function is defined such that it is zero at t0,
    ! rises to 1 smoothly, stay constant for tflat, and returns
    ! smoothly to zero.
    if(turn_drive_off) then
       epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
       S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
       S = S / (1-epsilon)
    else
       epsilon = 0.5*(tanh((t0-tL)/twL) + 1)
       S = 0.5*(tanh((t-tL)/twL) + 1) - epsilon
       S = S / (1-epsilon)
    endif
    if(S<0) then
       S = 0.
    endif
    return
  end subroutine PFenvelope


  subroutine setup_spline_nonunif_1D_periodic_aux( node_pos, N, buf, ibuf)
    sll_real64, dimension(:), pointer :: node_pos,buf
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), pointer :: a, cts
    sll_int32, dimension(:), pointer  :: ipiv,ibuf
    sll_int32 :: i
    a    => buf(1:3*N)
    cts  => buf(3*N+1:10*N)
    ipiv => ibuf(1:N)

    node_pos(-1)=node_pos(N-1)-(node_pos(N)-node_pos(0))
    node_pos(-2)=node_pos(N-2)-(node_pos(N)-node_pos(0))
    node_pos(N+1)=node_pos(1)+(node_pos(N)-node_pos(0))
    node_pos(N+2)=node_pos(2)+(node_pos(N)-node_pos(0))
    

    
    !fill a with mesh information
    do i=0,N-1
      !subdiagonal terms
      a(3*i+1)=(node_pos(i+1)-node_pos(i))*(node_pos(i+1)-node_pos(i))/((node_pos(i+1)-node_pos(i-1))*(node_pos(i+1)-node_pos(i-2)))
      !a(3*i+1)=(node_pos(i+1)-node_pos(i))/(node_pos(i+1)-node_pos(i-1))*((node_pos(i+1)-node_pos(i))/(node_pos(i+1)-node_pos(i-2)))
      !overdiagonal terms
      !a(3*i+3)=(node_pos(i)-node_pos(i-1))*(node_pos(i)-node_pos(i-1))/((node_pos(i+1)-node_pos(i-1))*(node_pos(i+2)-node_pos(i-1)))      
      a(3*i+3)=(node_pos(i)-node_pos(i-1))/(node_pos(i+1)-node_pos(i-1))*((node_pos(i)-node_pos(i-1))/(node_pos(i+2)-node_pos(i-1)))      
      !diagonal terms
      a(3*i+2)=1.0_f64-a(3*i+1)-a(3*i+3)
    enddo

    !initialize the tridiagonal solver
    call setup_cyclic_tridiag (a, N, cts, ipiv)
        
  end subroutine setup_spline_nonunif_1D_periodic_aux

  subroutine setup_spline_nonunif_1D_hermite_aux( node_pos, N, buf, ibuf)
    sll_real64, dimension(:), pointer :: node_pos,buf
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), pointer :: a, cts
    sll_int32, dimension(:), pointer  :: ipiv,ibuf
    sll_int32 :: i,Np
    Np=N+1
    a    => buf(1:3*Np)
    cts  => buf(3*Np+1:10*Np)
    ipiv => ibuf(1:Np)
    
    !symmetric case
    node_pos(-1)=2._f64*node_pos(0)-node_pos(1)
    node_pos(-2)=2._f64*node_pos(0)-node_pos(2)
    node_pos(N+1)=2.0_f64*node_pos(N)-node_pos(N-1)
    node_pos(N+2)=2.0_f64*node_pos(N)-node_pos(N-2)

    !triple point case
    node_pos(-1)=node_pos(0)
    node_pos(-2)=node_pos(0)
    node_pos(N+1)=node_pos(N)
    node_pos(N+2)=node_pos(N)

    
    !fill a with mesh information
    !a(1)=0.0_f64
    !a(2)=(node_pos(2)-node_pos(0))/(node_pos(1)+node_pos(2)-2._f64*node_pos(0))
    !a(3)=1.0_f64-a(2)
    a(1)=0.0_f64
    a(2)=(node_pos(2)-node_pos(0))/(node_pos(2)-node_pos(-1))
    a(3)=1.0_f64-a(2)
    do i=1,N-1
      !subdiagonal terms
      a(3*i+1)=(node_pos(i+1)-node_pos(i))**2/((node_pos(i+1)-node_pos(i-1))*(node_pos(i+1)-node_pos(i-2)))
      !overdiagonal terms
      a(3*i+3)=(node_pos(i)-node_pos(i-1))**2/((node_pos(i+1)-node_pos(i-1))*(node_pos(i+2)-node_pos(i-1)))      
      !diagonal terms
      a(3*i+2)=1.0_f64-a(3*i+1)-a(3*i+3)
    enddo
    !a(3*N+2)=(node_pos(N)-node_pos(N-2))/(2._f64*node_pos(N)-node_pos(N-1)-node_pos(N-2))
    !a(3*N+1)=1.0_f64-a(3*N+2)
    !a(3*N+3)=0.0_f64
    a(3*N+2)=(node_pos(N)-node_pos(N-2))/(node_pos(N+1)-node_pos(N-2))
    a(3*N+1)=1.0_f64-a(3*N+2)
    a(3*N+3)=0.0_f64

    !initialize the tridiagonal solver
    call setup_cyclic_tridiag (a, Np, cts, ipiv)
        
  end subroutine setup_spline_nonunif_1D_hermite_aux




  
  subroutine compute_spline_nonunif_1D_periodic_aux( f, N, buf, ibuf, coeffs )
    sll_real64, dimension(:), allocatable :: f
    sll_real64, dimension(:), pointer :: buf,coeffs
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), pointer :: cts!, a
    sll_int32, dimension(:), pointer  :: ipiv,ibuf
    !sll_real64 :: linf_err,tmp
    sll_int32 :: i
    !a    => buf(1:3*N) 
    cts  => buf(3*N+1:10*N)
    ipiv => ibuf(1:N)

    !compute the spline coefficients
    call solve_cyclic_tridiag_double( cts, ipiv, f, N, coeffs(0:N-1) )
    coeffs(-1) = coeffs(N-1)
    coeffs(N)  = coeffs(0)
    coeffs(N+1) = coeffs(1)
    
    !linf_err=0._f64
    !do i=1,N    
    !  tmp=a(3*(i-1)+1)*coeffs(i-2)+a(3*(i-1)+2)*coeffs(i-1)+a(3*(i-1)+3)*coeffs(i)-f(i)
    !  if(abs(tmp)>linf_err)then
    !    linf_err=abs(tmp)
    !  endif
    !enddo
    !print *,'error of compute_spline=',linf_err
  end subroutine compute_spline_nonunif_1D_periodic_aux

  subroutine compute_spline_nonunif_1D_hermite_aux( f, N, buf, ibuf, coeffs, lift )
    sll_real64, dimension(:), pointer :: f,buf,coeffs
    sll_int32, intent(in) :: N
    sll_real64, intent(in) :: lift(4,2)
    sll_real64, dimension(:), pointer :: cts
    sll_int32, dimension(:), pointer  :: ipiv,ibuf
    sll_int32 :: i,Np
    Np=N+1
    cts  => buf(3*Np+1:10*Np)
    ipiv => ibuf(1:Np)

    !compute the spline coefficients with inplace tridiagonal solver
    coeffs(1:N-1)=f(2:N)
    coeffs(0)=f(1)+lift(1,1)
    coeffs(N)=f(N+1)+lift(1,2)
    call solve_cyclic_tridiag_double( cts, ipiv, coeffs(0:N), Np, coeffs(0:N) )
    !coeffs(-1) = coeffs(1)+lift(2,1)
    !coeffs(N+1) = coeffs(N-1)+lift(2,2)
    
    coeffs(-1) = lift(3,1)*coeffs(0)+lift(4,1)*coeffs(1)+lift(2,1)
    coeffs(N+1) = lift(3,2)*coeffs(N-1)+lift(4,2)*coeffs(N)+lift(2,2)
    
    
  end subroutine compute_spline_nonunif_1D_hermite_aux



  subroutine interpolate_array_value_nonunif_aux( a_in, a_out, n, node_pos,coeffs,n_cells )
    sll_int32, intent(in) :: n,n_cells
    sll_real64, dimension(1:n), intent(in)  :: a_in
    !sll_real64, dimension(-1:n+1), intent(in)  :: node_pos
    sll_real64, dimension(1:n), intent(out) :: a_out
    sll_real64                         :: x
    !type(cubic_nonunif_spline_1D), pointer      :: spline
    sll_int32 :: i,j,shift=3
    sll_real64 ::xx
    sll_real64, dimension(:), pointer :: Xj
    sll_real64, dimension(:), pointer :: coef,coeffs,node_pos
    sll_real64 :: w(4)
    
    !do i=-1,n_cells+1
    !  print *,i,node_pos(i)
    !enddo
    
    xx=a_in(1)
    !xx = (x-spline%xmin)/(spline%xmax-spline%xmin)
    if(.not.((xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64 ))) then
      print *,'bad_value of x=',x!, 'xmin=', spline%xmin, 'xmax=', spline%xmax, xx,xx-1.0_f64
      print *,'in subroutine interpolate_array_value_nonunif()'
      STOP
    endif
    !SLL_ASSERT( (xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64 ))
    
    !localization of xx
    j=0
    if (xx==1.0_f64) then
      j = n_cells      
    else
      do while(node_pos(j).le.xx)
        j = j+1
      enddo
    endif
    
    do i=1,n
    
      x = a_in(i)
      if(.not.((x .ge. 0.0_f64) .and. (x .le. 1.0_f64 ))) then
        print *,'bad_value of a_in(',i,')=',a_in(i)!, 'xmin=', spline%xmin, 'xmax=', spline%xmax 
        print *,'in subroutine interpolate_array_value_nonunif()'
        STOP
      endif
      if (x==1.0_f64) then
        j = n_cells!spline%n_cells      
      else
        if(x.ge.xx) then
          do while(node_pos(j).le.x)
            j = j+1
          enddo
        else
          do while(node_pos(j).gt.x)
            j = j-1
          enddo
          j=j+1
        endif  
      endif
      xx=x
      Xj => node_pos(j-shift:)
      
      !print *,i,Xj(0+shift),xx,Xj(1+shift)
      !print *,i,Xj(-2+shift:3+shift)
      !stop
      if(.not.((xx .ge. Xj(0+shift)) .and. (xx .lt. Xj(1+shift)))) then
        if(xx.ne.1.0_f64) then
          print *,Xj(0+shift),xx,Xj(1+shift)
        
          stop
        endif  
      endif
      !SLL_ASSERT( (xx .ge. Xj(0+shift)) .and. (xx .le. Xj(1+shift) ))
      !SLL_ASSERT( (Xj(0+shift)==spline%node_positions(j-1)) .and. (Xj(1+shift)==spline%node_positions(j)))
    
      !compute weights
      w(1)=(Xj(shift+1)-xx)*(Xj(shift+1)-xx)*(Xj(shift+1)-xx)&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+1)-Xj(shift-2)));    
      w(2)=(Xj(shift+1)-xx)*(Xj(shift+1)-xx)*(xx-Xj(shift-2))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+1)-Xj(shift-2)));
      w(2)=w(2)+(Xj(shift+2)-xx)*(Xj(shift+1)-xx)*(xx-Xj(shift-1))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+2)-Xj(shift-1)));
      w(2)=w(2)+(Xj(shift+2)-xx)*(Xj(shift+2)-xx)*(xx-Xj(shift+0))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+2)-Xj(shift-1)));    
      w(3)=(Xj(shift+1)-xx)*(xx-Xj(shift-1))*(xx-Xj(shift-1))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+2)-Xj(shift-1)));
      w(3)=w(3)+(Xj(shift+2)-xx)*(xx-Xj(shift-1))*(xx-Xj(shift+0))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+2)-Xj(shift-1)));
      w(3)=w(3)+(Xj(shift+3)-xx)*(xx-Xj(shift+0))*(xx-Xj(shift+0))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+3)-Xj(shift+0)));    
      w(4)=(xx-Xj(shift+0))*(xx-Xj(shift+0))*(xx-Xj(shift+0))&
      &/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+3)-Xj(shift+0)));
    
      !coef => spline%coeffs(j-1:)
      !print *,i,xx,j,w(1),w(2),w(3),w(4)!,Xj(-2:3)
      !a_out(i) = w(0)*coef(0)+w(1)*coef(1)+w(2)*coef(2)+w(3)*coef(3)

      coef => coeffs(j-2:)
      a_out(i) = w(1)*coef(1)+w(2)*coef(2)+w(3)*coef(3)+w(4)*coef(4)
      !print *,i,xx,j,w(1),w(2),w(3),w(4),coef(1:4),a_out(i)
      !print *,Xj(shift-2:shift+3)
      !stop
    enddo
    !do i=1,n
    !  print *,i,a_in(i),a_out(i)
    !enddo  
  end subroutine interpolate_array_value_nonunif_aux

  subroutine compute_derivative_spl_per(f,N)
    !in f(1:N) out f(1:N+1)
    !we have to solve
    !df(i-1)+4df(i)+df(i+1)=3(f(i+1)-f(i-1)), i=1,N
    !for i=1, df(0)=df(N) for i=N, df(N+1)=df(1)
    sll_real64,dimension(:),allocatable::f,df,rhs
    sll_real64, dimension(:), pointer :: buf
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), pointer :: a, cts
    sll_int32, dimension(:), pointer  :: ipiv,ibuf
    sll_int32 :: i,im1,ip1
    sll_real64  ::err,dx,tmp
    
    allocate(buf(10*N))
    allocate(ibuf(N))    
    a    => buf(1:3*N)
    cts  => buf(3*N+1:10*N)
    ipiv => ibuf(1:N)
    allocate(df(1:N),rhs(1:N))

    do i=1,N
      im1=i-1
      ip1=i+1
      if(i==1)then
        im1=N
      endif
      if(i==N)then
        ip1=1
      endif
      rhs(i)=3._f64*(f(ip1)-f(im1))
    enddo
    

    
    !fill a with mesh information
    do i=0,N-1
      !subdiagonal terms
      a(3*i+1)=1._f64
      !overdiagonal terms
      a(3*i+3)=1._f64
      !diagonal terms
      a(3*i+2)=4._f64
    enddo

    !initialize the tridiagonal solver
    call setup_cyclic_tridiag (a, N, cts, ipiv)
    call solve_cyclic_tridiag_double( cts, ipiv, rhs(1:N), N, df(1:N) )

    !check that the result is ok
    err=0._f64
    do i=1,N
      im1=i-1
      ip1=i+1
      if(i==1)then
        im1=N
      endif
      if(i==N)then
        ip1=1
      endif
      tmp=df(im1)+4._f64*df(i)+df(ip1)-3._f64*(f(ip1)-f(im1))
      if(abs(tmp)>err)then
        err=abs(tmp)
      endif
    enddo
    if(err>1.e-12)then
      print *,'#bad inversion for compute_derivative_spl_per'
      print *,'#err=',err
      stop
    endif
    !print *,'#err for compute_derivative_spl_per',err
    f(1:N)=df(1:N)
    f(N+1)=f(1)
    deallocate(df,rhs,buf,ibuf)
  
  end subroutine compute_derivative_spl_per

  subroutine compute_derivative_spl_hrmt(f,N)
    !inout is f(-1:N+1): f(0:N)+2 derivatives (left: f(-1); right f(N+1))
    !we have to solve
    !df(i-1)+4df(i)+df(i+1)=3(f(i+1)-f(i-1)), i=1,N-1
    !for i=1, df(0)=f(-1) for i=N-1, df(N)=f(N+1)
    sll_real64,dimension(:),allocatable::f
    sll_real64,dimension(:),allocatable::df,rhs
    sll_real64, dimension(:), pointer :: buf
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), pointer :: a, cts
    sll_int32, dimension(:), pointer  :: ipiv,ibuf
    sll_int32 :: i,im1,ip1,Np
    sll_real64  ::err,dx,tmp
    
    Np=N-1  
    allocate(buf(10*Np))
    allocate(ibuf(Np))    
    a    => buf(1:3*Np)
    cts  => buf(3*Np+1:10*Np)
    ipiv => ibuf(1:Np)
    allocate(df(1:N-1),rhs(1:N-1))


    rhs(1)=3._f64*(f(2)-f(0))-f(-1)
    do i=2,N-2
      rhs(i)=3._f64*(f(i+1)-f(i-1))
    enddo
    rhs(N-1)=3._f64*(f(N)-f(N-2))-f(N+1)    
    
    !fill a with mesh information
    !a(2) a(3)   ...      a(1)
    !a(4) a(5) a(6)
    !     a(3i-2) a(3i-1) a(3i), i=2,N-2
    !a(3(N-1)) ... a(3(N-1)-2)  a(3(N-1)-1)
    a(1)=0._f64;a(2)=4._f64;a(3)=1._f64
    do i=2,N-2
      a(3*i-2)=1._f64;a(3*i-1)=4._f64;a(3*i)=1._f64
    enddo
    a(3*(N-1))=0._f64;a(3*(N-1)-2)=1._f64;a(3*(N-1)-1)=4._f64
    
    !initialize the tridiagonal solver
    call setup_cyclic_tridiag (a, N-1, cts, ipiv)
    call solve_cyclic_tridiag_double( cts, ipiv, rhs(1:N-1), N-1, df(1:N-1) )

    !check that the result is ok
    err=0._f64
    tmp=f(-1)+4._f64*df(1)+df(2)-3._f64*(f(2)-f(0))
    if(abs(tmp)>err)then
      err=abs(tmp)
    endif
    do i=2,N-2
      tmp=df(i-1)+4._f64*df(i)+df(i+1)-3._f64*(f(i+1)-f(i-1))
      if(abs(tmp)>err)then
        err=abs(tmp)
      endif
    enddo
    tmp=df(N-2)+4._f64*df(N-1)+f(N+1)-3._f64*(f(N)-f(N-2))
    if(abs(tmp)>err)then
      err=abs(tmp)
    endif

    if(err>1.e-13)then
      print *,'bad inversion for compute_derivative_spl_hrmt'
      stop
    endif
    !print *,'#err for compute_derivative_spl_hrmt',err
    f(1:N-1)=df(1:N-1)
    f(0)=f(-1)
    f(N)=f(N+1)
    deallocate(df,rhs,buf,ibuf)
  
  
  
  end subroutine compute_derivative_spl_hrmt







end module interp_non_unif_pp_class
