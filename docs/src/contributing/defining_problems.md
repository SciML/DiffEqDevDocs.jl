## Developing A New Problem

To develop a new problem, you need to make a new `DEProblem` and a new `DESolution`.
The `DEProblem` type should hold all of the mathematical information about the
problem, and the `DESolution` should hold all of the information for the solution.
Then all that is required is to define a `solve(::DEProblem,*Extra Mesh Things*;kwargs)`
which takes in the problem and returns a solution. To add plotting functionality,
add a plot recipe for the solution type to `/general/plotrecipes`. For testing
that the algorithm works, add a dispatch for `test_convergence` which makes
a `ConvergenceSimulation` type. This type already has a plot recipe, so
plotting functionality will already be embedded. This requires that your
problem can take in a true solution, and has a field `errors` which is a
dictionary of symbols for the different error estimates (L2,L infinity, etc.)
