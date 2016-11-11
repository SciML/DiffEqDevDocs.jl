# Contributor's Guide

So you're looking to help out DifferentialEquations.jl? We'd be happy to have
your help. It is recommended you first discuss with some of the developers
[on the Gitter channel](https://gitter.im/JuliaDiffEq/Lobby)
to make sure that you're up-to-date with current developments.

## Developing New Solver Algorithms

The easiest way to get started would be to add new solver algorithms. This is a
pretty simple task as there are tools which put you right into the "hot loop".
For example, take a look at the ODE solver code. The mode `solve(::ODEProblem,::AbstractArray)`
is glue code to a bunch of solver algorithms. The algorithms which are coded
in DifferentialEquations.jl can be found in ode_integrators.jl. For example,
take a look at the Midpoint method's implementation:

```julia
function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Midpoint,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
end
```

The available items are all unloaded from the `integrator` in the `@ode_preamble`.
`@ode_loopheader` and `@ode_loopfooter` macros are for exiting at max iterations,
and plugging into the Juno progressbar. These are all defined
using the `@def` macro (they essentially copy-paste the code from the line which
says `@def ode_loopheader begin ... end`). Note that the loopfooter code takes
care of the code for doing the adaptive timestepping. All that is required for
the adaptivity is that the algorithm computes an error estimate `EEst` each time,
save the value `utmp` to be what will replace `u` if the step is not rejected,
and add the algorithm's symbol is added to the dictionary `ODE_DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS`
in `ode_constants.jl`. If implicit solving is needed (via NLsolve),
add the algorithm's symbol to `DIFFERENTIALEQUATIONSJL_IMPLICITALGS` and the
conditional dependency will be supplied. Note that you may need more function
arguments. Use another method as a template.

When the solver is completed, add a call to the solver in the glue code
`solve(::ODEProblem,::AbstractArray)` (you will see all the others),
add the symbol for the algorithm to `DIFFERENTIALEQUATIONSJL_ALGORITHMS`, and
the order to `DIFFERENTIALEQUATIONSJL_ORDERS`. It's that quick! Lastly, add
your method to the convergence tests in the appropriate /test file.  Feel free
to implement any interesting or educational algorithm: they don't have to be
the fastest and it is always is useful to have such algorithms (like Simpson's method)
available for demonstration purposes.

Adding algorithms to the other problems is very similar.

### Extras

If the method is a FSAL method
then it needs to be set it in `DIFFERENTIALEQUATIONSJL_FSALALGS` and `fsalfirst`
should be defined before the loop, with `fsallast` what's pushed up to `fsalfirst`
upon a successful step. See `:DP5` for an example.

It's usually wise to dispatch onto Number separately since that uses `f(t,u)`
instead of `f(t,u,du)`. The dispatch is chosen by setting the `uType` and
`rateType`, usually to either `<:Number` or `<:AbstractArray` (though they
should be the same).

If tests fail due to units (i.e. SIUnits), don't worry. I would be willing to fix
that up. To do so, you have to make sure you keep separate your `rateType`s and
your `uType`s since the rates from `f` will have units of `u` but divided by
a unit of time. If you simply try to write these into `u`, the units part will
fail (normally you have to multiply by a ``dt``).
