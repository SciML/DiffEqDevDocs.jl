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
  `struct SSPRK22 <: OrdinaryDiffEqAlgorithm end`.
  Although it does not have the FSAL property, this is set to true since the derivative
  at the start and the end of the interval are used for the Hermite interpolation,
  and so this is FSAL'd so that way only a single extra function evaluation occurs
  over the whole integration. This is done in `alg_utils.jl` via
  `isfsal(alg::SSPRK22) = true`. Additionally, the order is set in the same
  file via `alg_order(alg::SSPRK22) = 2`.
- The algorithm `SSPRK22` is exported in `OrdinaryDiffEq.jl`.
- In `caches.jl`, the two cache types `SSPRK22Cache` (for inplace updates) and
  `SSPRK22ConstantCache` are defined, similarly to the other ones.
  Note: `u_cache(c::SSPRK22Cache) = ()` and
  `du_cache(c::SSPRK22Cache) = (c.k,c.du,c.fsalfirst)` return the parts of the
  modifiable cache that are changed if the size of the ODE changes.
- A new file `perform_step/ssprk_perform_step.jl` has been used for the new
  implementations. For both types of caches, the functions `initialize!`
  and `perform_step!` are defined there.
- Finally, tests are added. A new file `test/ode/ode_ssprk_tests.jl` is created
  and included in `tests/runtests.jl` via
  `@time @testset "SSPRK Tests" begin include("ode/ode_ssprk_tests.jl") end`.
- Additionally, regression tests for the dense output are added in
  `test/ode/ode_dense_tests.jl`.

For more details, refer to https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/pull/40

### Self-Contained Example

```julia
using OrdinaryDiffEq
import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm,OrdinaryDiffEqConstantCache,
      alg_order, alg_cache, initialize!, perform_step!, @muladd, @unpack,
      constvalue

struct Mead <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
export Mead
alg_order(alg::Mead) = 4

struct MeadConstantCache <: OrdinaryDiffEqConstantCache
  b1
  b2
  b3
  b4
  b5
  b6
  c2
  c3
  c4
  c5
  c6
end

function MeadConstantCache(T1,T2)
  b1 = T1(-0.15108370762927)
  b2 = T1(0.75384683913851)
  b3 = T1(-0.36016595357907)
  b4 = T1(0.52696773139913)
  b5 = T1(0)
  b6 = T1(0.23043509067071)
  c2 = T2(0.16791846623918)
  c3 = T2(0.48298439719700)
  c4 = T2(0.70546072965982)
  c5 = T2(0.09295870406537)
  c6 = T2(0.76210081248836)
  MeadConstantCache(b1,b2,b3,b4,b5,b6,c2,c3,c4,c5,c6)
end

alg_cache(alg::Mead,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = MeadConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

function initialize!(integrator, cache::MeadConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::MeadConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack b1,b2,b3,b4,b5,b6,c2,c3,c4,c5,c6 = cache
  k1=integrator.fsalfirst #f(uprev,p,t)
  k2=f(uprev+c2*dt*k1,p,t+c2*dt)
  k3=f(uprev+c3*dt*k2,p,t+c3*dt)
  k4=f(uprev+c4*dt*k3,p,t+c4*dt)
  k5=f(uprev+c5*dt*k4,p,t+c5*dt)
  k6=f(uprev+c6*dt*k5,p,t+c6*dt)
  u=uprev+dt*(b1*k1+b2*k2+b3*k3+b4*k4+b6*k6)
  integrator.fsallast = f(u,p,t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

f = ODEFunction((u,p,t)->1.01u,
            analytic = (u0,p,t) -> u0*exp(1.01t))
prob = ODEProblem(f,1.01,(0.0,1.0))
sol = solve(prob,Mead(),dt=0.1)

using Plots
plot(sol)
plot(sol,denseplot=false,plot_analytic=true)

using DiffEqDevTools
dts = (1/2) .^ (8:-1:1)
sim = test_convergence(dts,prob,Mead())
sim.𝒪est[:final]
plot(sim)

# Exanple of a good one!
sim = test_convergence(dts,prob,BS3())
sim.𝒪est[:final]
plot(sim)
```

### Adding new exponential algorithms

The exponential algorithms follow the same recipe as the general algorithms, but there
are automation utilities that make this easier. It is recommended that you refer to one
of the model algorithms for reference:

- For traditional exponential Runge-Kutta type methods (that come with a corresponding
  Butcher table), refer to `ETDRK2`.
- For adaptive exponential Rosenbrock type methods, refer to `Exprb32`.
- For exponential propagation iterative Runge-Kutta methods (EPIRK), refer to `EPIRK5P1`.

The first two classes support two modes of operation: operator caching and Krylov
approximation. The `perform_step!` method in `perform_step/exponential_rk_perform_step.jl`,
as a result, is split into two branches depending on whether `alg.krylov` is true. The
caching branch utilizes precomputed operators, which are calculated by the `expRK_operators`
method in `caches/linear_nonlinear_caches.jl`. Both `expRK_operators` and the `arnoldi`/`phiv`
methods in `perform_step!` comes from the
[ExponentialUtilities](https://github.com/JuliaDiffEq/ExponentialUtilities.jl) package.

The EPIRK methods can only use Krylov approximation, and unlike the previous two they use
the timestepping variant `phiv_timestep`. The timestepping method follows the convention
of Neisen & Wright, and can be toggled to use adaptation by `alg.adaptive_krylov`.

Although the exponential integrators (especially the in-place version) can seem complex, they
share similar structures. The infrastructure for the exisitng exponential methods utilize the
fact to reduce boilerplate code. In particular, the cache construction code in
`caches/linear_nonlinear_caches.jl` and the `initialize!` method in
`perform_step/exponential_rk_perform_step.jl` can be mostly automated and only `perform_step!`
needs implementing.

Finally, to construct tests for the new exponential algorithm, append the new algorithm to
the corresponding algorithm class in `test/linear_nonlinear_convergence_tests.jl` and
`test/linear_nonlinear_krylov_tests.jl`.
