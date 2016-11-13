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
