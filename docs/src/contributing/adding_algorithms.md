# Adding Algorithms

New algorithms can either be added by extending one of the current solver (or
add-on packages), or by contributing a new package to the organization. If
it's a new problem (a new PDE, a new type of differential equation, a new
subclass of problems for which special methods exist, etc.) then the problem
and solution types should be added to `DiffEqBase` first.

After the problem and solutions are defined, the `solve` method should be implemented.
It should take in keyword arguments which match the common interface (implement
"as many as possible"). One should note and document the amount of compatibility
with the common interface and Julia-defined types. After that, testing should be
done using `DiffEqDevTools`. Convergence tests and benchmarks should be included
to show the effectiveness of the algorithm and the correctness. Do not worry if
the algorithm is not "effective": the implementation can improve over time and
some algorithms useful just for the comparison they give!

After some development, one may want to document the algorithm in DiffEqBenchmarks
and DiffEqTutorials.


## Adding new algorithms to OrdinaryDiffEq

This recipe has been used to add the strong stability preserving Runge-Kutta methods
`SSPRK22`, `SSPRK33`, and `SSPRK104` to `OrdinaryDiffEq`. `SSPRK22` will be used
as an example.

- To create a new solver, two (three) types have to be created.
  The first is the algorithm `SSPRK22` used for dispatch, the other ones are
  the corresponding caches `SSPRK22Cache` (for inplace updates) and
  `SSPRK22ConstantCache`.
- The algorithm is defined in `algorithms.jl` as
  `immutable SSPRK22 <: OrdinaryDiffEqAlgorithm end`.
  Although it has not the FSAL property, this is set to true since the derivative
  at the start and the end of the interval are used for the Hermite interpolation,
  and so this is FSAL'd so that way only a single extra function evaluation occurs
  over the whole integration. This is done in `alg_utils.jl` via
  `isfsal(alg::SSPRK22) = true`. Additionally, the order is set in the same
  file via `alg_order(alg::SSPRK22) = 2`.
- The algorithm `SSPRK22`is exported in `OrdinaryDiffEq.jl`.
- In `caches.jl`, the two cache types `SSPRK22Cache` (for inplace updates) and
  `SSPRK22ConstantCache` are defined, similarly to the other ones.
  Note: `u_cache(c::SSPRK22Cache) = ()` and
  `du_cache(c::SSPRK22Cache) = (c.k,c.du,c.fsalfirst)` return the parts of the
  modifiable cache that are changed if the size of the ODE changes.
- A new file `integrators/ssprk_integrators.jl`has been used for the new
  implementations. For both types of caches, the functions `initialize!`
  and `perform_step!` are defined there.
- Finally, tests are added. A new file `test/ode/ode_ssprk_tests.jl` is created
  and included in `tests/runtests.jl` via 
  `@time @testset "SSPRK Tests" begin include("ode/ode_ssprk_tests.jl") end`.
- Additionally, regression tests for the dense output are added in
  `test/ode/ode_dense_tests.jl`.

For more details, refer to https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/pull/40
