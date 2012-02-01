PROGRAM POISSON3D_WITH__FFTW3

  use sll_collective


#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_remap.h"
  IMPLICIT NONE
  
  sll_int32 :: I, J, K, MX, MY, MZ
  sll_int32 :: NX, NY, NZ
  sll_int32 :: NCOUNT, T0, TN, ISTEP, NSTEP
  sll_int32 :: NTHREADS = 4
  sll_int32 :: FORWARD, BACKWARD
  sll_int64 :: col_size, myrank
  sll_int32 :: local_sz_X, local_sz_Y, local_sz_Z 
  sll_int32 :: NPX, NPY, NPZ ! Numbers of procs in the directions
  sll_int32 :: e, e1, e2, e3
  sll_int32 :: IERR, gI, gJ, gK

  sll_real64 :: X, Y, Z, TBEGIN, TEND, ERR
  sll_real64 :: DX, DY, DZ, TIME, DT, PI 
  sll_real64 :: CX, CY, CZ, VPX, VPY, VPZ 
  
  sll_real64, DIMENSION(:,:,:), ALLOCATABLE :: B, C, D, B_global
  sll_real64, DIMENSION(:,:,:), ALLOCATABLE :: U, F, G

  sll_int32, dimension(1:3) :: global

  type(layout_3D_t), pointer     :: layout1, layout2
  type(remap_plan_3D_t), pointer :: rmp3

  col_size = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  e = int(log(real(col_size))/log(2.))
  e1 = e/3
  e2 = (e-e1)/2
  e3 = e - (e1+e2)
  NPX = 2**e1
  NPY = 2**e2
  NPZ = 2**e3

  if ( col_size > min(NX-2,NY-2,NZ-2) ) then
     if (myrank==0) then
        print*, 'The number of processors must be <=', min(NX,NY,NZ)
     endif
     stop
  endif
  if ( .not.(is_power_of_two(col_size)) .or. (.not.is_power_of_two(int(NX,i64)-2)) .or. & 
       (.not.is_power_of_two(int(NY,i64)-2) ) .or. (.not.is_power_of_two(int(NZ,i64)-2)) ) then
     if (myrank==0) then
        print*, 'This need to run with 2 power numbers of processors and of points.'
     endif
     stop
  endif  
  
  NX = 130; NY = 130; NZ = 130
  DX = 1d0 / NX ; DY = 1d0 / NY ; DZ = 1d0 / NZ
  DT = 0.1; NSTEP = 10
  PI = 4d0 * DATAN(1d0)
  
  SLL_ALLOCATE(U(NX,NY,NZ), IERR); U = 0.
  SLL_ALLOCATE(F(NX,NY,NZ), IERR); F = 0.
  SLL_ALLOCATE(G(NX,NY,NZ), IERR); G = 0.
  
  CALL INIT_POISSON_FFTW(NX,NY,NZ)
  CALL SYSTEM_CLOCK(COUNT=T0, COUNT_RATE=NCOUNT)
  CALL CPU_TIME(TBEGIN)
  
  TIME = 0.0
  DO ISTEP = 1, NSTEP !Loop over time
     
     DO K = 1, NZ
        Z = (K-1)*DZ
        DO J = 1, NY
           Y = (J-1)*DY
           DO I = 1, NX
              X = (I-1)*DX
              F(I,J,K) = 12.*PI*PI*SIN(2*PI*X)*SIN(2*PI*Y)*SIN(2*PI*Z)*COS(2*PI*TIME)
              G(I,J,K) = SIN(2*PI*X)*SIN(2*PI*Y)*SIN(2*PI*Z)*COS(2*PI*TIME)
           ENDDO
        ENDDO
     ENDDO
     
     CALL SOLVE_POISSON_FFTW(U, F, NX, NY, NZ )
     
     WRITE(*,'("ISTEP = ",I4," TIME = ",G15.5," NX,NY,NZ =",3I4)') &
          ISTEP, TIME, NX, NY, NZ
     
     ERR = MAXVAL(U-G)
     PRINT  '(3X,"Norm(|ux - u_a|)               : ",2X,1PE10.3,//)', ERR
     
     TIME = TIME+DT
     
  END DO !Next time step
  
  CALL FINALIZE_POISSON_FFTW()
  
  CALL SYSTEM_CLOCK(COUNT=TN, COUNT_RATE=NCOUNT)
  WRITE(*,"(' ELAPSED TIME ', G15.5, ' s')") (TN - T0)/FLOAT(NCOUNT)
  CALL CPU_TIME(TEND)
  WRITE(*,"(' CPU TIME ', G15.5, ' s')") TEND-TBEGIN

  SLL_DEALLOCATE_ARRAY( U, IERR )
  SLL_DEALLOCATE_ARRAY( F, IERR )
  SLL_DEALLOCATE_ARRAY( G, IERR )
  
CONTAINS
  
  SUBROUTINE INIT_POISSON_FFTW(NX,NY,NZ)
#include "fftw3.f"
    sll_int32, INTENT(IN) :: NX,NY,NZ
    
    MX = NX-2; MY = NY-2; MZ = NZ-2
    
    CX = 1.0 / (DX*DX)
    CY = 1.0 / (DY*DY)
    CZ = 1.0 / (DZ*DZ)

    layout1  => new_layout_3D( sll_world_collective )  
    call initialize_layout_with_distributed_3D_array( MX, MY, MZ, NPX, NPY, NPZ, layout1 )  
    call compute_local_sizes( layout1, local_sz_X, local_sz_Y, local_sz_Z )   
    
    !RHS
    SLL_ALLOCATE(B(local_sz_X, local_sz_Y, local_sz_Z ), IERR)
    SLL_ALLOCATE(D(local_sz_x, local_sz_Y, local_sz_Z ), IERR)
    SLL_ALLOCATE(C(local_sz_X, local_sz_Y, local_sz_Z ), IERR)
    
    !Compute eigen values, build the matrix
    DO K=1,local_sz_Z
       DO J=1,local_sz_Y
          DO I=1,local_sz_X
             global = local_to_global_3D( layout1, (/I, J, K/))
             gI = global(1)
             gJ = global(2)
             gK = global(3)
             VPX=1.0-COS(FLOAT(gI)*PI/FLOAT(MX+1))
             VPY=1.0-COS(FLOAT(gJ)*PI/FLOAT(MY+1))
             VPZ=1.0-COS(FLOAT(gK)*PI/FLOAT(MZ+1))
             D(I,J,K)= 2.0 * (CZ*VPZ + CX*VPX + CY*VPY) &
                  * (8*(MX+1)*(MY+1)*(MZ+1))
          END DO
       END DO
    END DO
    
    !Initialize FFTs
    CALL DFFTW_INIT_THREADS(IERR)
    IF (IERR == 0) STOP 'FFTW CAN''T USE THREADS'
    CALL DFFTW_PLAN_WITH_NTHREADS(NTHREADS)
    CALL DFFTW_PLAN_R2R_3D(FORWARD,MX,MY,MZ,C,B, &
         FFTW_RODFT00, FFTW_RODFT00,FFTW_RODFT00,&
         FFTW_PATIENT )
    CALL DFFTW_PLAN_R2R_3D(BACKWARD,MX,MY,MZ,B,C, &
         FFTW_RODFT00, FFTW_RODFT00,FFTW_RODFT00,&
         FFTW_PATIENT )

    SLL_DEALLOCATE_ARRAY( B, IERR)
    SLL_DEALLOCATE_ARRAY( D, IERR)
    SLL_DEALLOCATE_ARRAY( C, IERR) 
    
  END SUBROUTINE INIT_POISSON_FFTW
  
  SUBROUTINE SOLVE_POISSON_FFTW(PSI, RHS, NX, NY, NZ)
#include "fftw3.f"
    sll_int32, INTENT(IN)  :: NX, NY, NZ
    REAL(8), INTENT(OUT) :: PSI(NX,NY,NZ)
    REAL(8), INTENT(IN)  :: RHS(NX,NY,NZ)
    
    !Boundary conditions
    DO K=1,local_sz_Z
       DO J=1,local_sz_Y
          DO I=1,local_sz_X
             global = local_to_global_3D( layout1, (/I, J, K/))
             gI = global(1)
             gJ = global(2)
             gK = global(3)
             C(I,J,K) = RHS(gI+1,gJ+1,gK+1)
             IF ( gI == 1) THEN
                C(I,J,K) = C(I,J,K) + CX * PSI( 1,gJ+1,gK+1)
             ELSE IF (I == MX) THEN
                C(I,J,K) = C(I,J,K) + CX * PSI(NX,gJ+1,gK+1)
             ELSE IF ( J == 1) THEN
                C(I,J,K) = C(I,J,K) + CY * PSI( gI+1,1,gK+1)
             ELSE IF (gJ == MY) THEN
                C(I,J,K) = C(I,J,K) + CY * PSI(gI+1,NY,gK+1)
             ELSE IF ( K == 1) THEN
                C(I,J,K) = C(I,J,K) + CZ * PSI(gI+1,gJ+1,1)
             ELSE IF (gK == MZ) THEN
                C(I,J,K) = C(I,J,K) + CZ * PSI(gI+1,gJ+1,NZ)
             END IF
          END DO
       END DO
    END DO
!Forward FFT on the right hand side term
    CALL DFFTW_EXECUTE_R2R(FORWARD,C,B)
    
    B = B / D 

    layout2  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( MX, MY, MZ, 1, 1, 1, layout2 )
    call get_local_sz( layout2, local_sz_X, local_sz_Y, local_sz_Z )  
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, B)
    SLL_ALLOCATE( B_global(MX, MY, MZ), IERR )
    call apply_remap_3D( rmp3, B, B_global )
    
    !Backward FFT on the solution
    CALL DFFTW_EXECUTE_R2R(BACKWARD,B_global,PSI(2:NX-1,2:NY-1,2:NZ-1))

    SLL_DEALLOCATE_ARRAY( B_global, IERR ) 
    
    RETURN
    
  END SUBROUTINE SOLVE_POISSON_FFTW
  
  SUBROUTINE FINALIZE_POISSON_FFTW()
#include "fftw3.f"
    
    CALL DFFTW_DESTROY_PLAN(FORWARD)
    CALL DFFTW_DESTROY_PLAN(BACKWARD)
    CALL DFFTW_CLEANUP_THREADS(IERR)
    
  END SUBROUTINE FINALIZE_POISSON_FFTW
 
  subroutine compute_local_sizes( layout, loc_sz_i, loc_sz_j, loc_sz_k )
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_3D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_3D_collective(layout))
    i_min = get_layout_3D_i_min( layout, my_rank )
    i_max = get_layout_3D_i_max( layout, my_rank )
    j_min = get_layout_3D_j_min( layout, my_rank )
    j_max = get_layout_3D_j_max( layout, my_rank )
    k_min = get_layout_3D_k_min( layout, my_rank )
    k_max = get_layout_3D_k_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
  end subroutine compute_local_sizes
  
END PROGRAM POISSON3D_WITH__FFTW3
