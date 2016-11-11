# The DiffEq Internals

The DiffEq solvers, OrdineryDiffEq, StochasticDiffEq, FiniteElementDiffEq, etc.
all follow a similar scheme which leads to rapid development and high performance.
This portion of the documentation explains how the algorithms are written.

## Developing New Solver Algorithms

The easiest way to get started would be to add new solver algorithms. This is a
pretty simple task as there are tools which put you right into the "hot loop".
For example, take a look at the ODE solver code. The mode `solve(::ODEProblem,::OrdinaryDiffEqAlgorithm)`
is glue code to a bunch of solver algorithms. The algorithms which are coded
in DifferentialEquations.jl can be found in ode_integrators.jl. For example,
take a look at the Midpoint method's implementation (without the function header):

```julia
  @ode_preamble
  halfdt::tType = dt/2
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      u = u + dt.*f(t+halfdt,u+halfdt.*f(t,u))
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
```

The available items are all unloaded from the `integrator` in the `@ode_preamble`.
`@ode_loopheader` and `@ode_loopfooter` macros are for exiting at max iterations,
and plugging into the Juno progressbar. These are all defined
using the `@def` macro (they essentially copy-paste the code from the line which
says `@def ode_loopheader begin ... end`). Note that the loopfooter code takes
care of the code for doing the adaptive timestepping. All that is required for
the adaptivity is that the algorithm computes an error estimate `EEst` each time,
save the value `utmp` to be what will replace `u` if the step is not rejected.
If implicit solving is needed (via NLsolve),
add the algorithm's symbol to `isimplicit` and the
conditional dependency will be supplied. Note that you may need more function
arguments. Use another method as a template.

It's that quick! Lastly, add your method to the convergence tests in the appropriate /test file.  
Feel free to implement any interesting or educational algorithm: they don't have to be
the fastest and it is always is useful to have such algorithms (like Simpson's method)
available for demonstration purposes.

Adding algorithms to the other problems is very similar.

## Extras

If the method is a FSAL method then it needs to be set via `isfsal` and `fsalfirst`
should be defined before the loop, with `fsallast` what's pushed up to `fsalfirst`
upon a successful step. See `:DP5` for an example.

It's usually wise to dispatch onto Number separately since that uses `f(t,u)`
instead of `f(t,u,du)`. The dispatch is chosen by setting the `uType` and
`rateType`, usually to either `<:Number` or `<:AbstractArray` (though they
should be the same).

If tests fail due to units (i.e. Unitful), don't worry. I would be willing to fix
that up. To do so, you have to make sure you keep separate your `rateType`s and
your `uType`s since the rates from `f` will have units of `u` but divided by
a unit of time. If you simply try to write these into `u`, the units part will
fail (normally you have to multiply by a ``dt``).
