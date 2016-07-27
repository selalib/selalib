UNIFORMLY ACCURATE FORWARD SEMI-LAGRANGIAN METHODS FOR HIGHLY OSCILLATORY VLASOV-POISSON EQUATION
==================================================================================================

This application is running the numerical simulation of a Vlasov-Poisson
equation modeling charged particles in a beam submitted to a highly
oscillatory external electric field. The numerical scheme is designed
to be uniformly accurate with respect to
the size of the fast time oscillations of the solution, which means
that no time step refinement is required to simulate the problem.

The scheme combines the forward semi-Lagrangian method with a class
of Uniformly Accurate (UA) time integrators to solve the characteristics.
These UA time integrators are derived by means of a two-scale
formulation of the characteristics, with the introduction of an
additional periodic variable.

You can run the numerical experiments changing the following parameters:

- `example` : example number described in the article.
- `eta_min`   : left boundary of the  mesh
- `eta_max`   : right boundary of the mesh
- `n`         : number of mesh points
- `k`         : time step
- `eps`       : ratio between the characteristic lengths in the transverse and the longitudinal directions
- `final_time` : simulation duration time
- `plot_enabled`: enable image outputs (must be false to reduce computation time)
- `ntau`: is the number of points for time discretization.
