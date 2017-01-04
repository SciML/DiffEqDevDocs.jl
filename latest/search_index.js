var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DiffEq-Developer-Documentation-1",
    "page": "Home",
    "title": "DiffEq Developer Documentation",
    "category": "section",
    "text": "This is the developer documentation for the DiffEq ecosystem. It explains the common interface and some the package internals to help developers contribute.If you have any questions, or just want to chat about solvers/using the package, please feel free to use the Gitter channel. For bug reports, feature requests, etc., please submit an issue. If you're interested in contributing, please see the Contributor's Guide."
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "The DiffEq ecosystem is built around the common interface. The common interface is a type-based interface where users define problems as a type, and solvers plug into the ecosystem by defining an algorithm to give a new dispatch tosolve(prob,alg;kwargs...)There is then an ecosystem of add-on components which use the common solver interface to add analysis tools for differential equations."
},

{
    "location": "index.html#Bleeding-Edge-1",
    "page": "Home",
    "title": "Bleeding Edge",
    "category": "section",
    "text": "This package suite is still under heavy development. If you are a power user and would like to try out the latest features, it is recommended you use the MetaDiffEq metapackage. To do so, use the following commands:Pkg.clone(\"https://github.com/tbreloff/MetaPkg.jl\") # Install MetaPkg\nusing MetaPkg\nmeta_add(\"MetaDiffEq\") # Adds all of the packages, even those unregistered\nmeta_checkout(\"MetaDiffEq\") # Checks out the master branch on all of the packagesNote that this is for power users who are familiar with Julia. If you are having issues, please contact Chris Rackauckas in  the Gitter channel."
},

{
    "location": "index.html#Contributing-to-the-Ecosystem-1",
    "page": "Home",
    "title": "Contributing to the Ecosystem",
    "category": "section",
    "text": "There are many ways to help the ecosystem. One way you can contribute is to give pull requests (PRs) to existing packages. Another way to contribute is to add your own package to the ecosystem. Adding your own package to the ecosystem allows you to keep executive control and licensing over your methods, but allows users of DifferentialEquations.jl to use your methods via the common interface, and makes your package compatible with the add-on tools (sensitivity analysis, parameter estimation, etc). Note that one is required to move their package to the JuliaDiffEq organization so that way common maintenance (such as fixing deprication warnings, updating tests to newer versions, and emergency fixes / disabling) can be allowed by JuliaDiffEq members. However, the lead developer of the package maintains administrative control, and thus any change to the core algorithms by other JuliaDiffEq members will only be given through PRs.Even if you don't have the time to contribute new solver algorithms or add-on tools, there's always ways to help! Improved plot recipes and new series recipes are always nice to add more default plots. It is always helpful to have benchmarks between different algorithms to see \"which is best\". Adding examples IJulia notebooks to DiffEqTutorials.jl is a good way to share knowledge about DifferentialEquations.jl. Also, please feel free to comb through the solvers and look for ways to make them more efficient. Lastly, the documentation could always use improvements. If you have any questions on how to help, just ask them in the Gitter!"
},

{
    "location": "index.html#Contributor-Guide-1",
    "page": "Home",
    "title": "Contributor Guide",
    "category": "section",
    "text": "Pages = [\n  \"contributing/ecosystem_overview.md\",\n  \"contributing/adding_algorithms.md\",\n  \"contributing/defining_problems.md\",\n  \"contributing/diffeq_internals.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Algorithm-Development-Tools-1",
    "page": "Home",
    "title": "Algorithm Development Tools",
    "category": "section",
    "text": "The following algorithm development tools are provided by DiffEqDevTools.jlPages = [\n  \"alg_dev/test_problems.md\",\n  \"alg_dev/convergence.md\",\n  \"alg_dev/benchmarks.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Internal-Documentation-1",
    "page": "Home",
    "title": "Internal Documentation",
    "category": "section",
    "text": "Pages = [\n  \"internals/fem_tools.md\",\n  \"internals/extras.md\",\n  \"internals/solver_helpers.md\",\n  \"internals/notes_on_algorithms.md\",\n  \"internals/tableaus.md\"\n]\nDepth = 2"
},

{
    "location": "contributing/ecosystem_overview.html#",
    "page": "Ecosystem Overview",
    "title": "Ecosystem Overview",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/ecosystem_overview.html#Ecosystem-Overview-1",
    "page": "Ecosystem Overview",
    "title": "Ecosystem Overview",
    "category": "section",
    "text": "So you're looking to help out DifferentialEquations.jl? We'd be happy to have your help. It is recommended you first discuss with some of the developers on the Gitter channel to make sure that you're up-to-date with current developments."
},

{
    "location": "contributing/ecosystem_overview.html#The-Common-Interface-1",
    "page": "Ecosystem Overview",
    "title": "The Common Interface",
    "category": "section",
    "text": "The DiffEq ecosystem is built around the common interface. This is the interface for the solvers:solve(prob,alg;kwargs...)and the standard methods for dealing with solutions. Ecosystem builds problem types for solvers to act on, and add-on components which use the solution types for higher-level analysis like parameter estimation and sensitivity analysis.One can add components at any of these levels to improve the functionality of the system as a whole."
},

{
    "location": "contributing/ecosystem_overview.html#Organizational-Setup-1",
    "page": "Ecosystem Overview",
    "title": "Organizational Setup",
    "category": "section",
    "text": "JuliaDiffEq is setup in a distributed manner to allow developers to retain authoritative control and licensing for their own packages/algorithms, yet contribute to the greater ecosystem. This gives a way for researchers to target a wide audience of users, but not have to fully contribute to public packages or be restricted in licensing. At the center of the ecosystem is DiffEqBase which holds the Problem, Solution, and Algorithm types (the algorithms are defined in DiffEqBase to be accessible by the default_algorithm function. One can opt out of this). Then there's the component solvers, which includes the *DiffEq packages (OrdinaryDiffEq, StochasticDiffEq, etc.) which implement different methods for solve. Then there are the add-on packages, such as the DiffEq* packages (DiffEqParamEstim, DiffEqDevTools) which add functionality to the Problem+solve setup. Lastly, there's DifferentialEquations.jl which is a metapackage which holds all of these pieces together as one cohesive unit.If one wants their package to officially join the ecosystem, it will need to be moved to the JuliaDiffEq organization so that maintenance can occur (but the core algorithms will only be managed by the package owners themselves). The Algorithm types can then be moved to DiffEqBase, and after testing the package will be added to the list packages exported by DifferentialEquations.jl and the corresponding documentation."
},

{
    "location": "contributing/adding_algorithms.html#",
    "page": "Adding Algorithms",
    "title": "Adding Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/adding_algorithms.html#Adding-Algorithms-1",
    "page": "Adding Algorithms",
    "title": "Adding Algorithms",
    "category": "section",
    "text": "New algorithms can either be added by extending one of the current solver (or add-on packages), or by contributing a new package to the organization. If it's a new problem (a new PDE, a new type of differential equation, a new subclass of problems for which special methods exist, etc.) then the problem and solution types should be added to DiffEqBase first.After the problem and solutions are defined, the solve method should be implemented. It should take in keyword arguments which match the common interface (implement \"as many as possible\"). One should note and document the amount of compatibility with the common interface and Julia-defined types. After that, testing should be done using DiffEqDevTools. Convergence tests and benchmarks should be included to show the effectiveness of the algorithm and the correctness. Do not worry if the algorithm is not \"effective\": the implementation can improve over time and some algorithms useful just for the comparison they give!After some development, one may want to document the algorithm in DiffEqBenchmarks and DiffEqTutorials."
},

{
    "location": "contributing/defining_problems.html#",
    "page": "Developing A New Problem",
    "title": "Developing A New Problem",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/defining_problems.html#Developing-A-New-Problem-1",
    "page": "Developing A New Problem",
    "title": "Developing A New Problem",
    "category": "section",
    "text": "New problems should be defined for new types of differential equations, new partial differential equations, and special subclasses of differential equations for which solvers can dispatch on for better performance (one such subclass are MassMatrixODEProblem which is a DAEProblem, but methods can specialize on the specific structure).To develop a new problem, you need to make a new DEProblem and a new DESolution. These types belong in DiffEqBase and should be exported. The DEProblem type should hold all of the mathematical information about the problem (including all of the meshing information in both space and time), and the DESolution should hold all of the information for the solution. Then all that is required is to define a solve(::DEProblem,alg;kwargs) which takes in the problem and returns a solution.If the problem makes sense in a heirarchy, one should define a promotion structure. For example, ODEProblems are DAEProblems, and so by defining this structure methods written for DAEs can solve ODEs by promotion (such as how DASSL provides a BDF method).To add plotting functionality, add a plot recipe for the solution type. For testing one should create a separate DETestProblem and DETestSolution which holds the analytical solution, and/or extend appxtrue! in DiffEqDevTools for error analysis. Then to check that the algorithm works, add a dispatch for test_convergence which makes a ConvergenceSimulation type. This type already has a plot recipe, so plotting functionality will already be embedded. This requires that your problem can take in a true solution, and has a field errors which is a dictionary of symbols for the different error estimates (L2,L infinity, etc.)After these steps, update the documentation to include the new problem types and the new associated solvers."
},

{
    "location": "contributing/diffeq_internals.html#",
    "page": "The DiffEq Internals",
    "title": "The DiffEq Internals",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/diffeq_internals.html#The-DiffEq-Internals-1",
    "page": "The DiffEq Internals",
    "title": "The DiffEq Internals",
    "category": "section",
    "text": "The DiffEq solvers, OrdineryDiffEq, StochasticDiffEq, FiniteElementDiffEq, etc. all follow a similar scheme which leads to rapid development and high performance. This portion of the documentation explains how the algorithms are written."
},

{
    "location": "contributing/diffeq_internals.html#Developing-New-Solver-Algorithms-1",
    "page": "The DiffEq Internals",
    "title": "Developing New Solver Algorithms",
    "category": "section",
    "text": "The easiest way to get started would be to add new solver algorithms. This is a pretty simple task as there are tools which put you right into the \"hot loop\". For example, take a look at the ODE solver code. The mode solve(::ODEProblem,::OrdinaryDiffEqAlgorithm) is glue code to a bunch of solver algorithms. The algorithms which are coded in DifferentialEquations.jl can be found in ode_integrators.jl. For example, take a look at the Midpoint method's implementation (without the function header):  @ode_preamble\n  halfdt::tType = dt/2\n  @inbounds for T in Ts\n    while t < T\n      @ode_loopheader\n      u = u + dt.*f(t+halfdt,u+halfdt.*f(t,u))\n      @ode_numberloopfooter\n    end\n  end\n  return u,t,timeseries,tsThe available items are all unloaded from the integrator in the @ode_preamble. @ode_loopheader and @ode_loopfooter macros are for exiting at max iterations, and plugging into the Juno progressbar. These are all defined using the @def macro (they essentially copy-paste the code from the line which says @def ode_loopheader begin ... end). Note that the loopfooter code takes care of the code for doing the adaptive timestepping. All that is required for the adaptivity is that the algorithm computes an error estimate EEst each time, save the value utmp to be what will replace u if the step is not rejected. If implicit solving is needed (via NLsolve), add the algorithm's symbol to isimplicit and the conditional dependency will be supplied. Note that you may need more function arguments. Use another method as a template.It's that quick! Lastly, add your method to the convergence tests in the appropriate /test file.   Feel free to implement any interesting or educational algorithm: they don't have to be the fastest and it is always is useful to have such algorithms (like Simpson's method) available for demonstration purposes.Adding algorithms to the other problems is very similar."
},

{
    "location": "contributing/diffeq_internals.html#Extras-1",
    "page": "The DiffEq Internals",
    "title": "Extras",
    "category": "section",
    "text": "If the method is a FSAL method then it needs to be set via isfsal and fsalfirst should be defined before the loop, with fsallast what's pushed up to fsalfirst upon a successful step. See :DP5 for an example.It's usually wise to dispatch onto Number separately since that uses f(t,u) instead of f(t,u,du). The dispatch is chosen by setting the uType and rateType, usually to either <:Number or <:AbstractArray (though they should be the same).If tests fail due to units (i.e. Unitful), don't worry. I would be willing to fix that up. To do so, you have to make sure you keep separate your rateTypes and your uTypes since the rates from f will have units of u but divided by a unit of time. If you simply try to write these into u, the units part will fail (normally you have to multiply by a dt)."
},

{
    "location": "alg_dev/test_problems.html#",
    "page": "Test Problems",
    "title": "Test Problems",
    "category": "page",
    "text": ""
},

{
    "location": "alg_dev/test_problems.html#Test-Problems-1",
    "page": "Test Problems",
    "title": "Test Problems",
    "category": "section",
    "text": "For every problem, there is an equivalent TestProblem which also has a field for the analytical solution, and a TestSolution which holds the analytical solution and calculates errors. This allows for easy testing/development, and works with the convregence simulation and benchmarking tools by default. If the solution was a TestProblem and thus has an analytical solution, we also havesol.u_analytic # timeseries of analytical solution\nsol.prob.analytic(t) # The analytic solution at time tavailable for further analysis."
},

{
    "location": "alg_dev/test_problems.html#No-Analytical-Solution-1",
    "page": "Test Problems",
    "title": "No Analytical Solution",
    "category": "section",
    "text": "However, in many cases the analytical solution cannot be found, and therefore one uses a low-tolerance calculation as a stand-in for a solution. The JuliaDiffEq ecosystem supports this through the TestSolution type in DiffEqDevTools. There are three constructors. The code is simple, so here it is:type TestSolution <: DESolution\n  t\n  u\n  interp\n  dense\nend\n(T::TestSolution)(t) = T.interp(t)\nTestSolution(t,u) = TestSolution(t,u,nothing,false)\nTestSolution(t,u,interp) = TestSolution(t,u,interp,true)\nTestSolution(interp::DESolution) = TestSolution(nothing,nothing,interp,true)This acts like a solution. When used in conjunction with apprxtrue:appxtrue(sol::AbstractODESolution,sol2::TestSolution)you can use it to build a TestSolution from a problem (like ODETestSolution) which holds the errors  If you only give it t and u, then it can only calculate the final error. If the TestSolution has an interpolation, it will define timeseries and dense errors.(Note: I would like it so that way the timeseries error will be calculated on the times of sol.t in sol2.t which would act nicely with tstops and when interpolations don't exist, but haven't gotten to it!)These can then be passed to other functionality. For example, the benchmarking functions allow one to set appxsol which is a TestSolution for the benchmark solution to calculate errors against, and error_estimate allows one to choose which error estimate to use in the benchmarking (defaults to :final)."
},

{
    "location": "alg_dev/test_problems.html#Related-Functions-1",
    "page": "Test Problems",
    "title": "Related Functions",
    "category": "section",
    "text": "DiffEqDevTools.appxtrue\nFiniteElementDiffEq.FEMSolutionTS"
},

{
    "location": "alg_dev/convergence.html#",
    "page": "Convergence Simulations",
    "title": "Convergence Simulations",
    "category": "page",
    "text": ""
},

{
    "location": "alg_dev/convergence.html#Convergence-Simulations-1",
    "page": "Convergence Simulations",
    "title": "Convergence Simulations",
    "category": "section",
    "text": "The convergence simulation type is useful for deriving order of convergence estimates from a group of simulations. This object will automatically assemble error vectors into a more useful manner and provide plotting functionality. Convergence estimates are also given by pair-wise estimates.One can automatically have DifferentialEquations.jl perform the error analysis by passing a ConvergenceSimulation a vector of solutions, or using one of the provided test_convergence functions. These will give order of convergence estimates and provide plotting functionality. This requires that the true solution was provided in the problem definition.ConvergenceSimulations can either be created by passing the constructor the appropriate solution array or by using one of the provided test_convergence functions."
},

{
    "location": "alg_dev/convergence.html#The-ConvergenceSimulation-Type-1",
    "page": "Convergence Simulations",
    "title": "The ConvergenceSimulation Type",
    "category": "section",
    "text": "A type which holds the data from a convergence simulation."
},

{
    "location": "alg_dev/convergence.html#Fields-1",
    "page": "Convergence Simulations",
    "title": "Fields",
    "category": "section",
    "text": "solutions::Array{<:DESolution}: Holds all the PdeSolutions.\nerrors: Dictionary of the error calculations. Can contain:\nh1Errors: Vector of the H1 errors.\nl2Errors: Vector of the L2 errors.\nmaxErrors: Vector of the nodal maximum errors.\nnode2Errors: Vector of the nodal l2 errors.\nN: The number of simulations.\nauxdata: Auxillary data of the convergence simluation. Entries can include:\ndts: The dt's in the simulations.\ndxs: The dx's in the simulations.\nμs: The CFL μ's in the simulations.\nνs: The CFL ν's in the simulations.\n𝒪est: Dictionary of order estimates. Can contain:\nConvEst_h1: The H1 error order of convergence estimate for the convergence simulation.  Generated via log2(error[i+1]/error[i]). Thus only valid if generated by halving/doubling  the dt/dx. If alternate scaling, modify by dividing of log(base,ConvEst_h1)\nConvEst_l2: The L2 error order of convergence estimate for the convergence simulation.  Generated via log2(error[i+1]/error[i]). Thus only valid if generated by halving/doubling  the dt/dx. If alternate scaling, modify by dividing of log(base,ConvEst_l2)\nConvEst_max: The nodal maximum error order of convergence estimate for the convergence simulation.  Generated via log2(error[i+1]/error[i]). Thus only valid if generated by halving/doubling  the dt/dx. If alternate scaling, modify by dividing of log(base,ConvEst_max)\nConvEst_node2: The nodal l2 error order of convergence estimate for the convergence simulation.  Generated via log2(error[i+1]/error[i]). Thus only valid if generated by halving/doubling  the dt/dx. If alternate scaling, modify by dividing of log(base,ConvEst_node2)\nconvergence_axis: The axis along which convergence is calculated. For example, if  we calculate the dt convergence, convergence_axis is the dts used in the calculation."
},

{
    "location": "alg_dev/convergence.html#Plot-Functions-1",
    "page": "Convergence Simulations",
    "title": "Plot Functions",
    "category": "section",
    "text": "The plot functionality is provided by a Plots.jl recipe. What is plotted is a line series for each calculated error along the convergence axis. To plot a convergence simulation, simply use:plot(sim::ConvergenceSimulation)All of the functionality (keyword arguments) provided by Plots.jl are able to be used in this command. Please see the Plots.jl documentation for more information."
},

{
    "location": "alg_dev/convergence.html#ODE-1",
    "page": "Convergence Simulations",
    "title": "ODE",
    "category": "section",
    "text": "test_convergence(dts::AbstractArray,prob::AbstractODEProblem)Tests the order of the time convergence of the given algorithm on the given problem solved over the given dts."
},

{
    "location": "alg_dev/convergence.html#Keyword-Arguments-1",
    "page": "Convergence Simulations",
    "title": "Keyword Arguments",
    "category": "section",
    "text": "T: The final time. Default is 1\nsave_timeseries: Denotes whether to save at every timeseries_steps steps. Default is true.\ntimeseries_steps: Denotes the steps to save at if save_timeseries=true. Default is 1\nalg: The algorithm to test.\ntableau: The tableau used for generic methods. Defaults to ODE_DEFAULT_TABLEAU."
},

{
    "location": "alg_dev/convergence.html#SDE-1",
    "page": "Convergence Simulations",
    "title": "SDE",
    "category": "section",
    "text": "test_convergence(dts::AbstractArray,prob::AbstractSDEProblem)Tests the strong order time convergence of the given algorithm on the given problem solved over the given dts."
},

{
    "location": "alg_dev/convergence.html#Keyword-Arguments-2",
    "page": "Convergence Simulations",
    "title": "Keyword Arguments",
    "category": "section",
    "text": "T: The final time. Default is 1\nnumMonte: The number of simulations for each dt. Default is 10000.\nsave_timeseries: Denotes whether to save at every timeseries_steps steps. Default is true.\ntimeseries_steps: Denotes the steps to save at if save_timeseries=true. Default is 1\nalg: The algorithm to test. Defaults to \"EM\"."
},

{
    "location": "alg_dev/convergence.html#Poisson-1",
    "page": "Convergence Simulations",
    "title": "Poisson",
    "category": "section",
    "text": "test_convergence(dxs::AbstractArray,prob::PoissonProblem)Tests the convergence of the solver algorithm on the given Poisson problem with dxs as given. Uses the square mesh [0,1]x[0,1]."
},

{
    "location": "alg_dev/convergence.html#Keyword-Arguments-3",
    "page": "Convergence Simulations",
    "title": "Keyword Arguments",
    "category": "section",
    "text": "solver: Which solver to use. Default is \"Direct\"."
},

{
    "location": "alg_dev/convergence.html#Heat-1",
    "page": "Convergence Simulations",
    "title": "Heat",
    "category": "section",
    "text": "test_convergence(dts::AbstractArray,dxs::AbstractArray,prob::AbstractHeatProblem,convergence_axis)Tests the convergence of the solver algorithm on the given Heat problem with the dts and dxs as given. Uses the square mesh [0,1]x[0,1]. The convergence axis is the axis along which convergence is calculated. For example, when testing dt convergence, convergence_axis = dts."
},

{
    "location": "alg_dev/convergence.html#Keyword-Arguments-4",
    "page": "Convergence Simulations",
    "title": "Keyword Arguments",
    "category": "section",
    "text": "T: The final time. Defaults to 1\nalg: The algorithm to test. Default is \"Euler\"."
},

{
    "location": "alg_dev/convergence.html#Utilities-1",
    "page": "Convergence Simulations",
    "title": "Utilities",
    "category": "section",
    "text": ""
},

{
    "location": "alg_dev/convergence.html#Order-Estimation-1",
    "page": "Convergence Simulations",
    "title": "Order Estimation",
    "category": "section",
    "text": "calc𝒪estimates(error::Vector{Number})`Computes the pairwise convergence estimate for a convergence test done by halving/doubling stepsizes vialog2(error[i+1]/error[i])Returns the mean of the convergence estimates."
},

{
    "location": "alg_dev/benchmarks.html#",
    "page": "Benchmark Suite",
    "title": "Benchmark Suite",
    "category": "page",
    "text": ""
},

{
    "location": "alg_dev/benchmarks.html#Benchmark-Suite-1",
    "page": "Benchmark Suite",
    "title": "Benchmark Suite",
    "category": "section",
    "text": "DiffernetialEquations.jl provides a benchmarking suite to be able to test the difference in error, speed, and efficiency between algorithms. DifferentialEquations.jl includes current benchmarking notebooks to help users understand the performance of the methods. These benchmarking notebooks use the included benchmarking suite. There are two parts to the benchmarking suite: shootouts and work-precision. The Shootout tests methods head-to-head for timing and error on the same problem. A WorkPrecision draws a work-precision diagram for the algorithms in question on the chosen problem."
},

{
    "location": "alg_dev/benchmarks.html#Using-the-Benchmarking-Notebooks-1",
    "page": "Benchmark Suite",
    "title": "Using the Benchmarking Notebooks",
    "category": "section",
    "text": "To use the benchmarking notebooks, IJulia is required. The commands are as follows:using IJulia\nnotebook(dir = Pkg.dir(\"DifferentialEquations\")*\"/benchmarks\")"
},

{
    "location": "alg_dev/benchmarks.html#Shootout-1",
    "page": "Benchmark Suite",
    "title": "Shootout",
    "category": "section",
    "text": "A shootout is where you compare between algorithms. For example, so see how different Runge-Kutta algorithms fair against each other, one can define a setup which is a dictionary of Symbols to Any, where the symbol is the keyword argument. Then you call ode_shootout on that setup. The code is as follows:tspan = [0,10]\nsetups = [Dict(:alg=>:DP5)\n          Dict(:abstol=>1e-3,:reltol=>1e-6,:alg=>:ode45) # Fix ODE to be normal\n          Dict(:alg=>:dopri5)]\nprob = DifferentialEquations.prob_ode_large2Dlinear\nnames = [\"DifferentialEquations\";\"ODE\";\"ODEInterface\"]\nshoot = ode_shootout(prob,tspan,setups;dt=1/2^(10),names=names)Note that keyword arguments applied to ode_shootout are applie dot every run, so in this example every run has the same starting timestep.  Here we explicitly chose names. If you don't, then the algorithm name is the default. This returns a Shootout type where which holds the times it took for each algorithm and the errors. Using these, it calculates the efficiency defnied as 1/(error*time), i.e. if the error is low or the run was quick then it's efficient. print(shoot) will show all of this information, and plot(shoot) will show the efficiencies of the algorithms in comparison to each other.For every benchmark function there is a special keyword numruns which controls the number of runs used in the time estimate. To be more precise, these functions by default run the algorithm 20 times on the problem and take the average time. This amount can be increased and decreased as needed.A ShootoutSet is a where you define a vector of probs and tspans and run a shootout on each of these values."
},

{
    "location": "alg_dev/benchmarks.html#WorkPrecision-1",
    "page": "Benchmark Suite",
    "title": "WorkPrecision",
    "category": "section",
    "text": "A WorkPrecision calculates the necessary componnets of a work-precision plot. This shows how time scales with the user chosen tolerances on a given problem. To make a WorkPrecision, you give it a vector of absolute and relative tolerances:abstols = 1./10.^(3:10)\nreltols = 1./10.^(3:10)\nwp = ode_workprecision(prob,tspan,abstols,reltols;alg=:DP5,name=\"Dormand-Prince 4/5\")If we want to plot many WorkPrecisions together in order to compare between algorithms, you can make a WorkPrecisionSet. To do so, you pass the setups into the function as well:wp_set = ode_workprecision_set(prob,tspan,abstols,reltols,setups;dt=1/2^4,numruns=2)\nsetups = [Dict(:alg=>:RK4);Dict(:alg=>:Euler);Dict(:alg=>:BS3);\n          Dict(:alg=>:Midpoint);Dict(:alg=>:BS5);Dict(:alg=>:DP5)]\nwp_set = ode_workprecision_set(prob,tspan,abstols,reltols,setups;dt=1/2^4,numruns=2)Both of these types have a plot recipe to produce a work-precision diagram, and a print which will show some relevant information."
},

{
    "location": "internals/fem_tools.html#",
    "page": "Internal Finite Element Tools",
    "title": "Internal Finite Element Tools",
    "category": "page",
    "text": ""
},

{
    "location": "internals/fem_tools.html#Internal-Finite-Element-Tools-1",
    "page": "Internal Finite Element Tools",
    "title": "Internal Finite Element Tools",
    "category": "section",
    "text": ""
},

{
    "location": "internals/fem_tools.html#Mesh-Tools-1",
    "page": "Internal Finite Element Tools",
    "title": "Mesh Tools",
    "category": "section",
    "text": "FiniteElementDiffEq.CFLν\nFiniteElementDiffEq.CFLμ"
},

{
    "location": "internals/fem_tools.html#Solver-Tools-1",
    "page": "Internal Finite Element Tools",
    "title": "Solver Tools",
    "category": "section",
    "text": "FiniteElementDiffEq.∇basis\nFiniteElementDiffEq.quadfbasis\nFiniteElementDiffEq.quadpts\nFiniteElementDiffEq.quadpts1\nFiniteElementDiffEq.assemblematrix\nFiniteElementDiffEq.∇u"
},

{
    "location": "internals/fem_tools.html#Error-Tools-1",
    "page": "Internal Finite Element Tools",
    "title": "Error Tools",
    "category": "section",
    "text": "FiniteElementDiffEq.getH1error\nFiniteElementDiffEq.getL2error"
},

{
    "location": "internals/extras.html#",
    "page": "Extra Functions",
    "title": "Extra Functions",
    "category": "page",
    "text": ""
},

{
    "location": "internals/extras.html#Extra-Functions-1",
    "page": "Extra Functions",
    "title": "Extra Functions",
    "category": "section",
    "text": "FiniteElementDiffEq.getNoise\nDiffEqBase.numparameters\nDiffEqBase.Tableau\nDiffEqBase.DEProblem"
},

{
    "location": "internals/solver_helpers.html#",
    "page": "Solver Extras",
    "title": "Solver Extras",
    "category": "page",
    "text": ""
},

{
    "location": "internals/solver_helpers.html#Solver-Extras-1",
    "page": "Solver Extras",
    "title": "Solver Extras",
    "category": "section",
    "text": ""
},

{
    "location": "internals/solver_helpers.html#SDE-Solver-Extras-1",
    "page": "Solver Extras",
    "title": "SDE Solver Extras",
    "category": "section",
    "text": "StochasticDiffEq.monteCarloSim\nStochasticDiffEq.RosslerSRI\nStochasticDiffEq.RosslerSRA\nStochasticDiffEq.constructSRA1\nStochasticDiffEq.constructSRIW1\nStochasticDiffEq.checkSRAOrder\nStochasticDiffEq.checkSRIOrder"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.GSδq!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.GSδq!",
    "category": "Function",
    "text": "GSδq!(δq,rp,Δxs)\n\nPerforms a Gauss-Seidel iteration for δq.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.GSu!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.GSu!",
    "category": "Function",
    "text": "GSu!(u,f₁,Δxs,p,ugD,grids,ux,uy)\n\nPerforms a Gauss-Seidel iteration on u.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.calc_rp!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.calc_rp!",
    "category": "Function",
    "text": "calc_rp!(rp,u,v,Δxs,g,px,py)\n\nCalculates the rp from the u and v's.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.update_p!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.update_p!",
    "category": "Function",
    "text": "update_p!(p,δq,Δxs)\n\nUpdates p given δq\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.update_v!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.update_v!",
    "category": "Function",
    "text": "update_v!(v,δq,Δxs)\n\nUpdates v given δq\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.uzawa_p!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.uzawa_p!",
    "category": "Function",
    "text": "uzawa_p!(p,u,v,Δxs,g,px,py)\n\nSolves for p from u and v using an Uzawa update.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.stokes_restriction",
    "page": "Solver Extras",
    "title": "StokesDiffEq.stokes_restriction",
    "category": "Function",
    "text": "stokes_restriction(u,v,p,Δxs,grids,mins,maxs,ugD,vgD)\n\nRestricts the Stokes problem to the coarsegrid.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.stokes_prolongation",
    "page": "Solver Extras",
    "title": "StokesDiffEq.stokes_prolongation",
    "category": "Function",
    "text": "stokes_prolongation(u,v,p,Δxs,grids,mins,maxs,ugD,vgD)\n\nProlongates the Stokes problem to the fine grid\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.update_u!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.update_u!",
    "category": "Function",
    "text": "update_u!(u,δq,Δxs)\n\nUpdates u given δq\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StokesDiffEq.GSv!",
    "page": "Solver Extras",
    "title": "StokesDiffEq.GSv!",
    "category": "Function",
    "text": "GSv!(v,f₂,Δxs,p,vgD,grids,vx,vy)\n\nPerforms a Gauss-Seidel iteration on v.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#Stationary-Stokes-1",
    "page": "Solver Extras",
    "title": "Stationary Stokes",
    "category": "section",
    "text": "StokesDiffEq.GSδq!\nStokesDiffEq.GSu!\nStokesDiffEq.calc_rp!\nStokesDiffEq.update_p!\nStokesDiffEq.update_v!\nStokesDiffEq.uzawa_p!\nStokesDiffEq.stokes_restriction\nStokesDiffEq.stokes_prolongation\nStokesDiffEq.update_u!\nStokesDiffEq.GSv!"
},

{
    "location": "internals/notes_on_algorithms.html#",
    "page": "Notes on Algorithms",
    "title": "Notes on Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "internals/notes_on_algorithms.html#Notes-on-Algorithms-1",
    "page": "Notes on Algorithms",
    "title": "Notes on Algorithms",
    "category": "section",
    "text": "This page is a supplemental page which details some facts about the chosen algorithms, why some I took the time to make optimized versions for, and for others why they were ignored."
},

{
    "location": "internals/notes_on_algorithms.html#Explicit-Runge-Kutta-ODE-Algorithms-1",
    "page": "Notes on Algorithms",
    "title": "Explicit Runge-Kutta ODE Algorithms",
    "category": "section",
    "text": "From what I can tell, this is by far the most comprehensive comparison of Explicit Runge-Kutta ODE algorithms that you'll find."
},

{
    "location": "internals/notes_on_algorithms.html#Implementations-1",
    "page": "Notes on Algorithms",
    "title": "Implementations",
    "category": "section",
    "text": "The different implementations have been benchmarked against each other. The efficiency was calculated by weighing both the time and error on classic test problems. To make clear distinctions, solver options were tweaked to many different settings, including:Matching errors\nMatching runtimes\nMatching settings\nLow/High toleranceThe DifferentialEquations.jl implementations of the explicit Runge-Kutta solvers are by a good margin the most efficient implementations of the given algorithms. They utilize many extra tricks, nice caching, and threading if available, to vastly outperform the other methods in terms of efficiency (even with threading disabled). :DP5 performs much better than :dopri5, which vastly outperform ode45 (whose stepsize algorithm tends to have issues on some quasi-stiff problems). :DP8 performs better than dop853 in some cases, worse in others. Both vastly outperform ode78.For this reason, the DifferentialEquations.jl non-stiff algorithms are the recommended implementations. ODEInterface non-stiff algorithms are only recommended for historical purposes (i.e. to match previous results) or to try dop853 on a problem (work is being to find out what the difference is and squash the competition here!). The ODE.jl algorithms are not recommended for any serious use (the package is essentially deprecated: it's slow, gets high error, the timestepping algorithm is not robust, and doesn't implement many methods)."
},

{
    "location": "internals/notes_on_algorithms.html#Order-4-1",
    "page": "Notes on Algorithms",
    "title": "Order 4-",
    "category": "section",
    "text": "At this stage, coefficient of the truncation error seems to win out, or you are willing to live with low tolerance anyways. Thus Bogacki-Shampine is the clear winner in this category because at order 2/3 with FASL it has minimal numbers of function evaluations but also is stable enough to step as needed. All other methods don't compare because of the FASL property boosting the order and thus the stability (for low orders, it pretty much holds that higher order = higher stability (for optimal number of steps), which is not true as we go higher), making it more stable and have less error for lower numbers of function evaluations than the others in this category."
},

{
    "location": "internals/notes_on_algorithms.html#Order-5-1",
    "page": "Notes on Algorithms",
    "title": "Order 5",
    "category": "section",
    "text": "[Note that for all of these Peter Stone's modifications do not seem to be helpful since, although they lower the truncation error, they also modify the stability region in ways that can be worrisome (mostly they shrink the stability in the complex axis near the origin, making the problems not as suitable for a \"general purpose default\" like one would hope with a 4/5 solver)]The \"clear choice\" is the Dormand-Prince 4/5 pair. This is the pair which is used by default as ode45 in MATLAB, and serves similar functions in scipy, ODE.jl, etc. The golden standard implementation is Hairer's DOPRI5 (offered by ODEInterface.jl). After optimizations, DifferentialEquations.jl's native DP5 solver is much more efficient (between 4x-400x) than DOPRI5's, with various design choices factoring into this (which are documented in the benchmarks). This is pre-threading, and within method threading will likely be at least doubled or tripled when threading is enabled. Thus it's clear that the reference implementation to try other methods against is the DifferentialEquations.jl DP5 method.It's obvious that anything before Dormand-Prince 4/5's pair is simply not as good because of the optimizations on the local truncation error coefficient and the fact that FASL schemes essentially have one less function evaluation. So the previous algorithms were implemented as tableaus for the historical reasons but dealt with no further. These methods include the Runge, Cassity, Butcher, Fehlburg, Lawson, Luther and Konen, and Kutta schemes.The next set of schemes are the Papakostas-Papageorgiou schemes. The problem is that they don't really get the much lower on the error than DP5, but also have wacky stability near the origin.Tsitouras's looks to be a good match against DP5 as a 6-stage scheme to take on DP5. Its stability is similar to DP5 but its first error term is an order of magnitude smaller. More tests will likely determine that this is much better than DP5 in accordance with his paper.Lastly, there are the 7-stage schemes. The more recent one is due to Sharp and Smart, but I am ignoring this because its error term is almost an order of magnitude larger than the BS pair, and its stability region is wonky near the origin. Its only plus over the BS pair is that it has a slightly larger stability in the real axis, which is not important when paired with adaptive stepping and for use on non-stiff problems.That leaves us with the Bogacki-Shampine pair. This pair gets more than an order of magnitude lower truncation error, enhanced complex stability, and two error estimators to make it more robust. In fact, this is the default which is chosen in Mathematica. Its downside is that since it is an 8-stage scheme, it requires an additional function evaluation.Further tests will likely narrow this down to Tsitouras vs Bogacki-Shampine. Who will come out on top? Who knows."
},

{
    "location": "internals/notes_on_algorithms.html#Order-6-1",
    "page": "Notes on Algorithms",
    "title": "Order 6",
    "category": "section",
    "text": "Sharp-Verner has bad complex stability near the origin. I don't like any of the Peter Stone modifications here. Butcher and Chummund methods have stability issues near the origin as well. Huta's method has too high of an error coefficient. Verner's 1991 has bad complex stability. Same as the most robust. The Verner \"most efficient\" has really good stability and error coefficient. In fact, nothing is even close except for Tsitouras' method. The DP method is two orders of magnitude higher in error coefficient than Verner. The Luther methods have too much error. Same as Tsitouras-Papakostas and  M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki.Without a doubt the winner is the Verner \"most efficient\"."
},

{
    "location": "internals/notes_on_algorithms.html#Order-7-1",
    "page": "Notes on Algorithms",
    "title": "Order 7",
    "category": "section",
    "text": "The Enright-Verner and other Verner methods all have stability issues near the origin in the complex plane and higher error coefficients. Sharp and Smart have higher error coefficients. Peter Stone's methods all have higher error. It's very clear that the best here is the Tanaka-Yamashita (efficient, not the stable) method by far."
},

{
    "location": "internals/notes_on_algorithms.html#Order-8-1",
    "page": "Notes on Algorithms",
    "title": "Order 8",
    "category": "section",
    "text": "The Cooper-Verner methods do not have an error estimate and so no adaptive timestepping can be done. This is a deal-breaker. Going into this one would think that the clear winner would be Dormand-Prince 8. But that's not the case. However, that's comparing the classical 1981 DP87. Notice that the code for Dop853 is based off of the 1989 paper which has different coefficients (and currently I have no analysis for this).The other methods include Verner's Maple dverk78 which is bested in both stability and error coefficient by Enright-Verner which is bested by Tsitouras-Papakostas.Thus the final showdown is between DP853 vs the Tsitouras-Papakostas pair."
},

{
    "location": "internals/notes_on_algorithms.html#Order-9-1",
    "page": "Notes on Algorithms",
    "title": "Order 9",
    "category": "section",
    "text": "The Tsitouras scheme and the Sharp scheme have funky stability near the origin. Verner's schemes are much safer, and with similar error. They clearly dominate this category."
},

{
    "location": "internals/notes_on_algorithms.html#Order-10-1",
    "page": "Notes on Algorithms",
    "title": "Order 10",
    "category": "section",
    "text": "Curtis' scheme has more function evaluations than needed, and Peter Stone's modification reduces the truncation error by a lot but adds three more function evaluations. Thus Hairer's 17 stage scheme (whose error and stability is similar to Curtis') is clearly better. Once again Peter Stone's modification adds three steps but does not reduce the truncation error here, so the unmodified version does better.Tom Baker's method increases the stability region to something which is more than necessary but adds 4 function evaluations to do so (without lowering the error very much). Ono's scheme minimizes the error more than Hairer's here, with all else being basically the same. The Peter Stone methods add a lot of function evaluations (5+) and so they would only be useful in the case where the function evaluations are quick yet you still want really small error. Even then I'm not convinced they are better than the other methods, or better than the higher order methods which use less steps. The stability is only okay.The Feagin scheme is fine, but with more error and less stability than the Hairer scheme. Thus it seems clear that Hairer's method dominates this category. However, that's only because it does not include an error estimate. Feagin's scheme is close in error and stability, but includes an error estimate which can be used for adaptivity, making it the choice in this category."
},

{
    "location": "internals/notes_on_algorithms.html#Order-11-1",
    "page": "Notes on Algorithms",
    "title": "Order 11",
    "category": "section",
    "text": "The order 11 schemes are due to Tom Baker at the University of Teeside. They have a nice sparsity pattern and receive slightly lower truncation error coefficents than the Feagin, but Feagin's dominates by being \"almost order 13\" anyways so while a nice try the order 11 scheme is likely overwhelmed in any case where it would be useful."
},

{
    "location": "internals/notes_on_algorithms.html#Order-12-1",
    "page": "Notes on Algorithms",
    "title": "Order 12",
    "category": "section",
    "text": "Here there are the Feagin schemes and Ono's scheme. Ono's scheme gets horrible stability with more error and so it's not in the running. Peter Stone's modifications do not make a substantive change, and where they do they get rid of the nice property that the Feagin 12 method satisfies many of the higher order conditions as well, making it look even higher order on some problems. Thus the standard Feagin 12 seems to win out in this category."
},

{
    "location": "internals/notes_on_algorithms.html#Order-14-1",
    "page": "Notes on Algorithms",
    "title": "Order 14",
    "category": "section",
    "text": "In this category there is just the Feagin. Peter Stone's modification barely changes anything in the analysis so I did not even attempt it."
},

{
    "location": "internals/tableaus.html#",
    "page": "ODE Tableaus",
    "title": "ODE Tableaus",
    "category": "page",
    "text": ""
},

{
    "location": "internals/tableaus.html#ODE-Tableaus-1",
    "page": "ODE Tableaus",
    "title": "ODE Tableaus",
    "category": "section",
    "text": ""
},

{
    "location": "internals/tableaus.html#Explicit-Runge-Kutta-Methods-1",
    "page": "ODE Tableaus",
    "title": "Explicit Runge-Kutta Methods",
    "category": "section",
    "text": "constructEuler - Euler's 1st order method.\nconstructHuen() Huen's order 2 method.\nconstructRalston() - Ralston's order 2 method.\nconstructKutta3 - Kutta's classic 3rd order method\nconstructRK4 - The classic 4th order \"Runge-Kutta\" method\nconstructRK438Rule - The classic 4th order \"3/8th's Rule\" method\nconstructBogakiShampine3() - Bogakai-Shampine's 2/3 method.\nconstructRKF4() - Runge-Kutta-Fehlberg 3/4.\nconstructRKF5() - Runge-Kutta-Fehlberg 4/5.\nconstructRungeFirst5() - Runge's first 5th order method.\nconstructCassity5() - Cassity's 5th order method.\nconstructLawson5() - Lawson's 5th order method.\nconstructLutherKonen5 - Luther-Konen's first 5th order method.\nconstructLutherKonen52() - Luther-Konen's second 5th order method.\nconstructLutherKonen53() - Luther-Konen's third 5th order method.\nconstructPapakostasPapaGeorgiou5() - Papakostas and PapaGeorgiou more stable order 5 method.\nconstructPapakostasPapaGeorgiou52() - Papakostas and PapaGeorgiou more efficient order 5 method.\nconstructTsitouras5() - Tsitouras's order 5 method.\nconstructBogakiShampine5() - Bogaki and Shampine's Order 5 method.\nconstructSharpSmart5() - Sharp and Smart's Order 5 method.\nconstructCashKarp() - Cash-Karp method 4/5.\nconstructDormandPrince() - Dormand-Prince 4/5.\nconstructButcher6() - Butcher's first order 6 method.\nconstructButcher62() - Butcher's second order 6 method.\nconstructButcher63() - Butcher's third order 6 method.\nconstructDormandPrince6() - Dormand-Prince's 5/6 method.\nconstructSharpVerner6() Sharp-Verner's 5/6 method.\nconstructVerner916() - Verner's more efficient order 6 method (1991).\nconstructVerner9162() - Verner's second more efficient order 6 method (1991).\nconstructVernerRobust6() - Verner's \"most robust\" order 6 method.\nconstructVernerEfficient6() - Verner's \"most efficient\" order 6 method.\nconstructPapakostas6() - Papakostas's order 6 method.\nconstructLawson6() - Lawson's order 6 method.\nconstructTsitourasPapakostas6() - Tsitouras and Papakostas's order 6 method.\nconstructDormandLockyerMcCorriganPrince6() - the Dormand-Lockyer-McCorrigan-Prince order 6 method.\nconstructTanakaKasugaYamashitaYazaki6A() - Tanaka-Kasuga-Yamashita-Yazaki order 6 method A.\nconstructTanakaKasugaYamashitaYazaki6B() - Tanaka-Kasuga-Yamashita-Yazaki order 6 method B.\nconstructTanakaKasugaYamashitaYazaki6C() - Tanaka-Kasuga-Yamashita-Yazaki order 6 method C.\nconstructTanakaKasugaYamashitaYazaki6D() - Tanaka-Kasuga-Yamashita-Yazaki order 6 method D.\nconstructMikkawyEisa() - Mikkawy and Eisa's order 6 method.\nconstructChummund6() - Chummund's first order 6 method.\nconstructChummund62() - Chummund's second order 6 method.\nconstructHuta6() - Huta's first order 6 method.\nconstructHuta62() - Huta's second order 6 method.\nconstructVerner6() - An old order 6 method attributed to Verner.\nconstructDverk() - The classic DVERK algorithm attributed to Verner.\nconstructClassicVerner6() - A classic Verner order 6 algorithm (1978).\nconstructButcher7() - Butcher's order 7 algorithm.\nconstructClassicVerner7()- A classic Verner order 7 algorithm (1978).\nconstructVernerRobust7() - Verner's \"most robust\" order 7 algorithm.\nconstructTanakaYamashitaStable7() - Tanaka-Yamashita more stable order 7 algorithm.\nconstructTanakaYamashitaEfficient7() - Tanaka-Yamashita more efficient order 7 algorithm.\nconstructSharpSmart7() - Sharp-Smart's order 7 algorithm.\nconstructSharpVerner7() - Sharp-Verner's order 7 algorithm.\nconstructVerner7() - Verner's \"most efficient\" order 7 algorithm.\nconstructVernerEfficient7() - Verner's \"most efficient\" order 7 algorithm.\nconstructClassicVerner8() - A classic Verner order 8 algorithm (1978).\nconstructCooperVerner8() - Cooper-Verner's first order 8 algorithm.\nconstructCooperVerner82() - Cooper-Verner's second order 8 algorithm.\nconstructTsitourasPapakostas8() - Tsitouras-Papakostas order 8 algorithm.\nconstructdverk78() - The classic order 8 DVERK algorithm.\nconstructEnrightVerner8() - Enright-Verner order 8 algorithm.\nconstructCurtis8() - Curtis' order 8 algorithm.\nconstructVerner8() - Verner's \"most efficient\" order 8 algorithm.\nconstructRKF8() - Runge-Kutta-Fehlberg Order 7/8 method.\nconstructDormandPrice8() - Dormand-Prince Order 7/8 method.\nconstructDormandPrince8_64bit() - Dormand-Prince Order 7/8 method. Coefficients are rational approximations good for 64 bits.\nconstructVernerRobust9() - Verner's \"most robust\" order 9 method.\nconstructVernerEfficient9() - Verner's \"most efficient\" order 9 method.\nconstructSharp9() - Sharp's order 9 method.\nconstructTsitouras9() - Tsitouras's first order 9 method.\nconstructTsitouras92() - Tsitouras's second order 9 method.\nconstructCurtis10() - Curtis' order 10 method.\nconstructOno10() - Ono's order 10 method.\nconstructFeagin10Tableau() - Feagin's order 10 method.\nconstructCurtis10() - Curtis' order 10 method.\nconstructBaker10() - Baker's order 10 method.\nconstructHairer10() Hairer's order 10 method.\nconstructFeagin12Tableau() - Feagin's order 12 method.\nconstructOno12() - Ono's order 12 method.\nconstructFeagin14Tableau() Feagin's order 14 method."
},

{
    "location": "internals/tableaus.html#Implicit-Runge-Kutta-Methods-1",
    "page": "ODE Tableaus",
    "title": "Implicit Runge-Kutta Methods",
    "category": "section",
    "text": "constructImplicitEuler - The 1st order Implicit Euler method.\nconstructMidpointRule - The 2nd order Midpoint method.\nconstructTrapezoidalRule - The 2nd order Trapezoidal rule (2nd order LobattoIIIA)\nconstructLobattoIIIA4 - The 4th order LobattoIIIA\nconstructLobattoIIIB2 - The 2nd order LobattoIIIB\nconstructLobattoIIIB4 - The 4th order LobattoIIIB\nconstructLobattoIIIC2 - The 2nd order LobattoIIIC\nconstructLobattoIIIC4 - The 4th order LobattoIIIC\nconstructLobattoIIICStar2 - The 2nd order LobattoIIIC*\nconstructLobattoIIICStar4 - The 4th order LobattoIIIC*\nconstructLobattoIIID2 - The 2nd order LobattoIIID\nconstructLobattoIIID4 - The 4th order LobattoIIID\nconstructRadauIA3 - The 3rd order RadauIA\nconstructRadauIA5 - The 5th order RadauIA\nconstructRadauIIA3 - The 3rd order RadauIIA\nconstructRadauIIA5 - The 5th order RadauIIA"
},

{
    "location": "internals/tableaus.html#Base.length-Tuple{DiffEqBase.ODERKTableau}",
    "page": "ODE Tableaus",
    "title": "Base.length",
    "category": "Method",
    "text": "Base.length(tab::ODERKTableau)\n\nDefines the length of a Runge-Kutta method to be the number of stages.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.stability_region",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.stability_region",
    "category": "Function",
    "text": "stability_region(z,tab::ODERKTableau)\n\nCalculates the stability function from the tableau at z. Stable if <1.\n\nr(z) = fracdet(I-zA+zeb^T)det(I-zA)\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqBase.ExplicitRKTableau",
    "page": "ODE Tableaus",
    "title": "DiffEqBase.ExplicitRKTableau",
    "category": "Type",
    "text": "ExplicitRKTableau\n\nHolds a tableau which defines an explicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqBase.ImplicitRKTableau",
    "page": "ODE Tableaus",
    "title": "DiffEqBase.ImplicitRKTableau",
    "category": "Type",
    "text": "ImplicitRKTableau\n\nHolds a tableau which defines an implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqBase.ODERKTableau",
    "page": "ODE Tableaus",
    "title": "DiffEqBase.ODERKTableau",
    "category": "Type",
    "text": "ODERKTableau: A Runge-Kutta Tableau for an ODE integrator\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.ODE_DEFAULT_TABLEAU",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.ODE_DEFAULT_TABLEAU",
    "category": "Constant",
    "text": "ODE_DEFAULT_TABLEAU\n\nSets the default tableau for the ODE solver. Currently Dormand-Prince 4/5.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#Tableau-Methods-1",
    "page": "ODE Tableaus",
    "title": "Tableau Methods",
    "category": "section",
    "text": "Base.length(::DiffEqBase.ODERKTableau)\nDiffEqDevTools.stability_region\nDiffEqBase.ExplicitRKTableau\nDiffEqBase.ImplicitRKTableau\nDiffEqBase.ODERKTableau\nOrdinaryDiffEq.ODE_DEFAULT_TABLEAU"
},

{
    "location": "internals/tableaus.html#Explicit-Tableaus-1",
    "page": "ODE Tableaus",
    "title": "Explicit Tableaus",
    "category": "section",
    "text": "DiffEqDevTools.constructEuler\nDiffEqDevTools.constructRalston\nDiffEqDevTools.constructHeun\nDiffEqDevTools.constructKutta3\nOrdinaryDiffEq.constructBS3\nDiffEqDevTools.constructBogakiShampine3\nDiffEqDevTools.constructRK4\nDiffEqDevTools.constructRK438Rule\nDiffEqDevTools.constructRKF4\nDiffEqDevTools.constructRKF5\nDiffEqDevTools.constructCashKarp\nDiffEqDevTools.constructDormandPrince\nOrdinaryDiffEq.constructBS5\nDiffEqDevTools.constructPapakostasPapaGeorgiou5\nDiffEqDevTools.constructPapakostasPapaGeorgiou52\nDiffEqDevTools.constructTsitouras5\nDiffEqDevTools.constructLutherKonen5\nDiffEqDevTools.constructLutherKonen52\nDiffEqDevTools.constructLutherKonen53\nDiffEqDevTools.constructRungeFirst5\nDiffEqDevTools.constructLawson5\nDiffEqDevTools.constructSharpSmart5\nDiffEqDevTools.constructBogakiShampine5\nDiffEqDevTools.constructCassity5\nDiffEqDevTools.constructButcher6\nDiffEqDevTools.constructButcher62\nDiffEqDevTools.constructButcher63\nDiffEqDevTools.constructVernerRobust6\nDiffEqDevTools.constructTanakaKasugaYamashitaYazaki6A\nDiffEqDevTools.constructTanakaKasugaYamashitaYazaki6B\nDiffEqDevTools.constructTanakaKasugaYamashitaYazaki6C\nDiffEqDevTools.constructTanakaKasugaYamashitaYazaki6D\nDiffEqDevTools.constructHuta6\nDiffEqDevTools.constructHuta62\nDiffEqDevTools.constructVerner6\nDiffEqDevTools.constructDormandPrince6\nDiffEqDevTools.constructSharpVerner6\nDiffEqDevTools.constructVern6\nDiffEqDevTools.constructClassicVerner6\nDiffEqDevTools.constructChummund6\nDiffEqDevTools.constructChummund62\nDiffEqDevTools.constructPapakostas6\nDiffEqDevTools.constructLawson6\nDiffEqDevTools.constructTsitourasPapakostas6\nDiffEqDevTools.constructDormandLockyerMcCorriganPrince6\nDiffEqDevTools.constructVernerEfficient6\nDiffEqDevTools.constructMikkawyEisa\nDiffEqDevTools.constructVernerEfficient7\nDiffEqDevTools.constructClassicVerner7\nDiffEqDevTools.constructSharpVerner7\nDiffEqDevTools.constructTanakaYamashitaStable7\nDiffEqDevTools.constructSharpSmart7\nDiffEqDevTools.constructTanakaYamashitaEfficient7\nDiffEqDevTools.constructVernerRobust7\nOrdinaryDiffEq.constructTanYam7\nDiffEqDevTools.constructEnrightVerner7\nDiffEqDevTools.constructDormandPrince8\nDiffEqDevTools.constructRKF8\nDiffEqDevTools.constructCooperVerner8\nDiffEqDevTools.constructCooperVerner82\nDiffEqDevTools.constructTsitourasPapakostas8\nDiffEqDevTools.constructEnrightVerner8\nDiffEqDevTools.constructdverk78\nDiffEqDevTools.constructClassicVerner8\nDiffEqDevTools.constructDormandPrince8_64bit\nDiffEqDevTools.constructCurtis8\nOrdinaryDiffEq.constructTsitPap8\nDiffEqDevTools.constructSharp9\nDiffEqDevTools.constructTsitouras9\nDiffEqDevTools.constructTsitouras92\nDiffEqDevTools.constructVernerEfficient9\nOrdinaryDiffEq.constructVern9\nDiffEqDevTools.constructVerner916\nDiffEqDevTools.constructVerner9162\nDiffEqDevTools.constructVernerRobust9\nDiffEqDevTools.constructFeagin10\nDiffEqDevTools.constructFeagin10Tableau\nDiffEqDevTools.constructOno10\nDiffEqDevTools.constructCurtis10\nDiffEqDevTools.constructHairer10\nDiffEqDevTools.constructBaker10\nDiffEqDevTools.constructFeagin12\nDiffEqDevTools.constructOno12\nDiffEqDevTools.constructFeagin12Tableau\nDiffEqDevTools.constructFeagin14\nDiffEqDevTools.constructFeagin14Tableau"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructImplicitEuler",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructImplicitEuler",
    "category": "Function",
    "text": "Implicit Euler Method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructMidpointRule",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructMidpointRule",
    "category": "Function",
    "text": "Order 2 Midpoint Method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTrapezoidalRule",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTrapezoidalRule",
    "category": "Function",
    "text": "Order 2 Trapezoidal Rule (LobattoIIIA2)\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIIA4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIIA4",
    "category": "Function",
    "text": "LobattoIIIA Order 4 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIIB2",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIIB2",
    "category": "Function",
    "text": "LobattoIIIB Order 2 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIIB4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIIB4",
    "category": "Function",
    "text": "LobattoIIIB Order 4 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIIC2",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIIC2",
    "category": "Function",
    "text": "LobattoIIIC Order 2 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIIC4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIIC4",
    "category": "Function",
    "text": "LobattoIIIC Order 4 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIICStar2",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIICStar2",
    "category": "Function",
    "text": "LobattoIIIC* Order 2 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIICStar4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIICStar4",
    "category": "Function",
    "text": "LobattoIIIC* Order 4 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIID2",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIID2",
    "category": "Function",
    "text": "LobattoIIID Order 2 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLobattoIIID4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLobattoIIID4",
    "category": "Function",
    "text": "LobattoIIID Order 4 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructGL2",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructGL2",
    "category": "Function",
    "text": "Gauss-Legendre Order 2.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructGL4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructGL4",
    "category": "Function",
    "text": "Gauss-Legendre Order 4.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructGL6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructGL6",
    "category": "Function",
    "text": "Gauss-Legendre Order 6.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRadauIA3",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRadauIA3",
    "category": "Function",
    "text": "RadauIA Order 3 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRadauIA5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRadauIA5",
    "category": "Function",
    "text": "RadauIA Order 5 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRadauIIA3",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRadauIIA3",
    "category": "Function",
    "text": "RadauIIA Order 3 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRadauIIA5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRadauIIA5",
    "category": "Function",
    "text": "RadauIIA Order 5 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#Implicit-Tableaus-1",
    "page": "ODE Tableaus",
    "title": "Implicit Tableaus",
    "category": "section",
    "text": "DiffEqDevTools.constructImplicitEuler\nDiffEqDevTools.constructMidpointRule\nDiffEqDevTools.constructTrapezoidalRule\nDiffEqDevTools.constructLobattoIIIA4\nDiffEqDevTools.constructLobattoIIIB2\nDiffEqDevTools.constructLobattoIIIB4\nDiffEqDevTools.constructLobattoIIIC2\nDiffEqDevTools.constructLobattoIIIC4\nDiffEqDevTools.constructLobattoIIICStar2\nDiffEqDevTools.constructLobattoIIICStar4\nDiffEqDevTools.constructLobattoIIID2\nDiffEqDevTools.constructLobattoIIID4\nDiffEqDevTools.constructGL2\nDiffEqDevTools.constructGL4\nDiffEqDevTools.constructGL6\nDiffEqDevTools.constructRadauIA3\nDiffEqDevTools.constructRadauIA5\nDiffEqDevTools.constructRadauIIA3\nDiffEqDevTools.constructRadauIIA5"
},

]}
