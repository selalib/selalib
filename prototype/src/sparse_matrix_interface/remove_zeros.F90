subroutine remove_zeros ( NR, NC, A, JA, IA, ANNZ)
! A, JA, IA must be pointers
implicit none
        INTEGER :: NR
        INTEGER :: NC
        INTEGER :: ANNZ        
        REAL(8), DIMENSION(ANNZ)  :: A
        INTEGER, DIMENSION(ANNZ) :: JA
        INTEGER, DIMENSION(NR+1) :: IA
        ! LOCAL VARIABLES
        REAL(8), DIMENSION(:), POINTER  :: T
        INTEGER, DIMENSION(:), POINTER :: IT
        INTEGER, DIMENSION(:), POINTER :: JT
        integer  :: TNNZ
        integer  :: li_i
        integer  :: li_k
        integer  :: li_j
        REAL(8) :: epsilon        

        epsilon = 1.0e-13

        ALLOCATE(T(ANNZ))
        ALLOCATE(JT(ANNZ))
        ALLOCATE(IT(NR+1))

        !... parcours pour copier les elts non nuls dans le nouveau tableau
        !... nbr d elt non nuls
        TNNZ   = 0                
        do li_i = 1 , NR
        
                IT(li_i) = TNNZ + 1
                
                do li_k = IA ( li_i ) , IA ( li_i + 1 ) - 1
                
                        li_j = JA ( li_k )
                        
                        if ( dabs ( A ( li_k ) ) >=  epsilon ) then
                        
                                TNNZ = TNNZ + 1
                                
                                JT ( TNNZ ) = li_j
                                
                                T ( TNNZ )  =  A ( li_k )
                                
                        end if                        
                        
                end do
                
        end do        

        IT ( NR + 1 ) = TNNZ + 1                                        

        A = 0.0
        JA = 0
        IA = 0

        !... copy data from the temp matrix
        ANNZ = TNNZ
        IA ( 1 : NR + 1 ) = IT ( 1 : NR + 1 )
        JA ( 1 : TNNZ    ) = JT ( 1 : TNNZ )
        A  ( 1 : TNNZ    ) = T  ( 1 : TNNZ )        

        DEALLOCATE(T, IT, JT)
end subroutine remove_zeros
