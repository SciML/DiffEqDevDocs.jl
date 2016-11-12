# Developing A New Problem

To develop a new problem, you need to make a new `DEProblem` and a new `DESolution`.
These types belong in DiffEqBase and should be exported.
The `DEProblem` type should hold all of the mathematical information about the
problem (including all of the meshing information in both space and time),
and the `DESolution` should hold all of the information for the solution.
Then all that is required is to define a `solve(::DEProblem,alg;kwargs)`
which takes in the problem and returns a solution. To add plotting functionality,
add a plot recipe for the solution type. For testing one should create a
separate `DETestProblem` and `DETestSolution` which holds the analytical
solution, and/or extend `appxtrue!` in DiffEqDevTools for error analysis.
Then to check that the algorithm works, add a dispatch for `test_convergence`
which makes a `ConvergenceSimulation` type. This type already has a plot recipe, so
plotting functionality will already be embedded. This requires that your
problem can take in a true solution, and has a field `errors` which is a
dictionary of symbols for the different error estimates (L2,L infinity, etc.)

After these steps, update the documentation to include the new problem types and
the new associated solvers.
