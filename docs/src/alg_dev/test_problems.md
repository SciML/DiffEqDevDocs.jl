If the solution was a `TestProblem` and thus has an analytical solution, we also have

```julia
sol.u_analytic # timeseries of analytical solution
sol.prob.analytic(t) # The analytic solution at time t
```
