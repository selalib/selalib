.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _poisson_2d_cartesian:

2D Poisson solver in Cartesian coordinates
==========================================

This model implements the 2D poisson problem, in cartesian coordinates system, on a domain :math:`\Omega` with homogeous Dirichlet boundary conditions, for a given function :math:`f`.

.. math::

  - \nabla^2 u = f, ~ ~ ~ \Omega
    \\
    u = 0, ~~~ \partial \Omega

The variational formulation writes

.. math::

   \sum_j U^j \int_{\Omega} V_R^i V_R^j+  V_Z^i V_Z^j  = \int_{\Omega} f V_0^i 

which leads to the following linear system

.. math::

   \mathcal{S} U = F

where

.. math::
 
   \mathcal{S}_{ij} = \int_{\Omega} V_R^i V_R^j ~~~ \mbox{and} ~~~ F_i = \int_{\Omega} f V_0^i  

Test Case
_________

Only one single testcase is implemented, for a square domain. 

The analytical solution is

.. math::

   u(R,Z) = \sin(2 \pi m_1 R) \sin(2 \pi n_1 Z)

which leads to the following right hands 

.. math::

   f = - \nabla^2 u = 4 \pi^2 (m_1^2 + n_1^2) \sin(2 \pi m_1 R) \sin(2 \pi n_1 Z)

where :math:`m_1` and :math:`n_1` can be changed in the *parameter* input file using the variables **INT_MODES_M1** and **INT_MODES_N1**. 

Numerical results for the convergence order are given below

.. include:: ../../tests/test_1/results.rst

.. Local Variables:
.. mode: rst
.. End:
