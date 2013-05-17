!> @file sll_integration.F90
!! @namespace sll_integration
!! @author Madaule Eric
!! @brief Gauss-Lobatto and Gauss-Legendre integration tools

!! @details Here are several of the Gauss-Lobatto tools :\\
!!            ·Gauss-Lobatto points and weight,\\
!!            ·Gauss-Lobatto bases functions and the integral of their product,\\
!!            ·integral of product of Gauss-Lobatto function and their derivative.
!!
!!          The mass matrix (which is the integral of \phi_i \times \phi_j) is simply 
!!          diag(weigh), so there is no need to store it more than just the weigh.
!!
!!          We also need the derivative matrix D.
!!          \f[ D_{i,j}=\int \pih_i \phi^'_j f]
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and so have time to write it).
!!         
!> There is a low-level mathematical utility that applies the 
!> Gauss-Legendre method to compute numeric integrals.
!> This module aims at providing a single interface to the process of 
!> integrating a function on a given interval.
