# The Parameters Interface

## Building a New Problem with New Parameters

From a `DEProblem`, the function

```julia
problem_new_parameters(prob::DEProblem,p)
```

can be used to build a new problem type which uses the parameter vector `p`.
For problems which have parameters in multiple places, such as a `SDEProblem`
or a `MonteCarloProblem`, this `p` is split between these locations using
`num_params` to find out how much of `p` goes to each location.
