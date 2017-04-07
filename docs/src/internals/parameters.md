# The Parameters Interface

The parameters interface allows for parameters in functions to be added and
accessed by the various routines which needs parameters, and allows them
to be "buried" in the routines which need to not have explicit parameters.

## The Functions on Parameterized Functions

```julia
param_values(f) # An array of the parameters in the function
num_params(f) # The number of parameters in the function
```

## Building a New Problem with New Parameters

From a `DEProblem`, the function

```julia
problem_new_parameters(prob::DEProblem,p)
```

can be used to build a new problem type which uses the parameter vector `p`.
For problems which have parameters in multiple places, such as a `SDEProblem`
or a `MonteCarloProblem`, this `p` is split between these locations using
`num_params` to find out how much of `p` goes to each location.
