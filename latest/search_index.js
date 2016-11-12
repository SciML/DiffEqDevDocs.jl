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
    "location": "index.html#Other-Help-1",
    "page": "Home",
    "title": "Other Help",
    "category": "section",
    "text": "Even if you don't have the time to contribute new solver algorithms, there's always ways to help! Improved plot recipes and new series recipes are always nice to add more default plots. It is always helpful to have benchmarks between different algorithms to see \"which is best\". Adding examples IJulia notebooks to DiffEqTutorials.jl is a good way to share knowledge about DifferentialEquations.jl. Also, please feel free to comb through the solvers and look for ways to make them more efficient. Lastly, the documentation could always use improvements. If you have any questions on how to help, just ask them in the Gitter!"
},

{
    "location": "index.html#Contributor-Guide-1",
    "page": "Home",
    "title": "Contributor Guide",
    "category": "section",
    "text": "Pages = [\n  \"contributing/contributors_guide.md\",\n  \"contributing/adding_algorithms.md\",\n  \"contributing/defining_problems.md\",\n]\nDepth = 2"
},

{
    "location": "index.html#Algorithm-Development-Tools-1",
    "page": "Home",
    "title": "Algorithm Development Tools",
    "category": "section",
    "text": "The following algorithm development tools are provided by DiffEqDevTools.jlPages = [\n  \"alg_dev/convergence.md\",\n  \"alg_dev/benchmarks.md\"\n]\nDepth = 2"
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
    "page": "Contributor's Guide",
    "title": "Contributor's Guide",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/ecosystem_overview.html#Contributor's-Guide-1",
    "page": "Contributor's Guide",
    "title": "Contributor's Guide",
    "category": "section",
    "text": "So you're looking to help out DifferentialEquations.jl? We'd be happy to have your help. It is recommended you first discuss with some of the developers on the Gitter channel to make sure that you're up-to-date with current developments."
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
    "text": ""
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
    "text": "To develop a new problem, you need to make a new DEProblem and a new DESolution. The DEProblem type should hold all of the mathematical information about the problem, and the DESolution should hold all of the information for the solution. Then all that is required is to define a solve(::DEProblem,*Extra Mesh Things*;kwargs) which takes in the problem and returns a solution. To add plotting functionality, add a plot recipe for the solution type to /general/plotrecipes. For testing that the algorithm works, add a dispatch for test_convergence which makes a ConvergenceSimulation type. This type already has a plot recipe, so plotting functionality will already be embedded. This requires that your problem can take in a true solution, and has a field errors which is a dictionary of symbols for the different error estimates (L2,L infinity, etc.)"
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
    "location": "internals/fem_tools.html#FiniteElementDiffEq.CFLν",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.CFLν",
    "category": "Function",
    "text": "CFLν(dt,dx)`\n\nComputes the CFL-condition = dtdx\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.CFLμ",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.CFLμ",
    "category": "Function",
    "text": "CFLμ(dt,dx)`\n\nComputes the CFL-condition = dt(dx*dx)\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#Mesh-Tools-1",
    "page": "Internal Finite Element Tools",
    "title": "Mesh Tools",
    "category": "section",
    "text": "FiniteElementDiffEq.CFLν\nFiniteElementDiffEq.CFLμ"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.∇basis",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.∇basis",
    "category": "Function",
    "text": "∇basis(node,elem)\n\nReturns u of the barycentric basis elements.\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.quadfbasis",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.quadfbasis",
    "category": "Function",
    "text": "quadfbasis(f,gD,gN,A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars;gNquad𝒪=2)\n\nPerforms the order 2 quadrature to calculate the vector from the term fv for linear elements.\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.quadpts",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.quadpts",
    "category": "Function",
    "text": "quadpts(𝒪)\n\nReturns the quadrature points and ω's for and 𝒪  in 2D.\n\nReference: David Dunavant. High degree efficient symmetrical Gaussian quadrature rules for the triangle. International journal for numerical methods in engineering. 21(6):1129–1148, 1985.\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.quadpts1",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.quadpts1",
    "category": "Function",
    "text": "quadpts1(𝒪)\n\nReferences: Pavel Holoborodko: http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.assemblematrix",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.assemblematrix",
    "category": "Function",
    "text": "assemblematrix(node,elem;lumpflag=false,K=[])\n\nAssembles the stiffness matrix A as an approximation to Δ on the finite element mesh (node,elem). Also generates the mass matrix M. If lumpflag=true, then the mass matrix is lumped resulting in a diagonal mass matrix. Specify a diffusion constant along the nodes via K.\n\nReturns\n\nA = Stiffness Matrix\nM = Mass Matrix\narea = A vector of the calculated areas for each element.\n\n\n\nassemblematrix(FEMmesh::FEMmesh;lumpflag=false,K=[])\n\nAssembles the stiffness matrix A as an approximation to Δ on the finite element mesh (node,elem). Also generates the mass matrix M. If lumpflag=true, then the mass matrix is lumped resulting in a diagonal mass matrix. Specify a diffusion constant along the nodes via K.\n\nReturns\n\nA = Stiffness Matrix\nM = Mass Matrix\narea = A vector of the calculated areas for each element.\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.∇u",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.∇u",
    "category": "Function",
    "text": "∇u(node,elem,u,Dλ=[])\n\nEstimates u on the mesh (node,elem)\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#Solver-Tools-1",
    "page": "Internal Finite Element Tools",
    "title": "Solver Tools",
    "category": "section",
    "text": "FiniteElementDiffEq.∇basis\nFiniteElementDiffEq.quadfbasis\nFiniteElementDiffEq.quadpts\nFiniteElementDiffEq.quadpts1\nFiniteElementDiffEq.assemblematrix\nFiniteElementDiffEq.∇u"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.getH1error",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.getH1error",
    "category": "Function",
    "text": "function getH1error(node,elem,Du,uh,K=[],quad𝒪=[])\n\ngetH1error(fem_mesh::FEMmesh,Du,u)\n\nEstimates the H1 error between uexact and uh on the mesh (node,elem). It reads the mesh to estimate the element type and uses this to choose a quadrature 𝒪 unless specified. If K is specified then it is the diffusion coefficient matrix.\n\n\n\n"
},

{
    "location": "internals/fem_tools.html#FiniteElementDiffEq.getL2error",
    "page": "Internal Finite Element Tools",
    "title": "FiniteElementDiffEq.getL2error",
    "category": "Function",
    "text": "getL2error(node,elem,uexact,uh,quad𝒪=[])\n\ngetL2error(fem_mesh::FEMmesh,sol,u)\n\nEstimates the L2 error between uexact and uh on the mesh (node,elem). It reads the mesh to estimate the element type and uses this to choose a quadrature 𝒪 unless specified.\n\n\n\n"
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
    "location": "internals/extras.html#FiniteElementDiffEq.getNoise",
    "page": "Extra Functions",
    "title": "FiniteElementDiffEq.getNoise",
    "category": "Function",
    "text": "getNoise(N,node,elem;noisetype=:White)\n\nReturns a random vector corresponding to the noise type which was chosen.\n\n\n\n"
},

{
    "location": "internals/extras.html#DiffEqBase.numparameters",
    "page": "Extra Functions",
    "title": "DiffEqBase.numparameters",
    "category": "Function",
    "text": "numparameters(f)\n\nReturns the number of parameters of f for the method which has the most parameters.\n\n\n\n"
},

{
    "location": "internals/extras.html#DiffEqBase.Tableau",
    "page": "Extra Functions",
    "title": "DiffEqBase.Tableau",
    "category": "Type",
    "text": "Tableau: Holds the information for a Runge-Kutta Tableau\n\n\n\n"
},

{
    "location": "internals/extras.html#DiffEqBase.DEProblem",
    "page": "Extra Functions",
    "title": "DiffEqBase.DEProblem",
    "category": "Type",
    "text": "DEProblem: Defines differential equation problems via its internal functions\n\n\n\n"
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
    "location": "internals/solver_helpers.html#StochasticDiffEq.monteCarloSim",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.monteCarloSim",
    "category": "Function",
    "text": "monteCarloSim(dt::Number,prob::SDEProblem)\n\nPerforms a parallel Monte-Carlo simulation to solve the SDE problem with dt numMonte times. Returns a vector of solution objects.\n\nKeyword Arguments\n\nT - Final time. Default is 1.\nnumMonte - Number of Monte-Carlo simulations to run. Default is 10000\nsave_timeseries - Denotes whether save_timeseries should be turned on in each run. Default is false.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StochasticDiffEq.RosslerSRI",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.RosslerSRI",
    "category": "Type",
    "text": "RosslerSRI\n\nHolds the Butcher tableaus for a Rosser SRI method.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StochasticDiffEq.RosslerSRA",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.RosslerSRA",
    "category": "Type",
    "text": "RosslerSRA\n\nHolds the Butcher tableaus for a Rosser SRA method.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StochasticDiffEq.constructSRA1",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.constructSRA1",
    "category": "Function",
    "text": "constructSRA1()\n\nConstructs the taleau type for the SRA1 method.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StochasticDiffEq.constructSRIW1",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.constructSRIW1",
    "category": "Function",
    "text": "constructSRIW1()\n\nConstructs the tableau type for the SRIW1 method.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StochasticDiffEq.checkSRAOrder",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.checkSRAOrder",
    "category": "Function",
    "text": "checkSRAOrder(RosslerSRI)\n\nDetermines whether the order conditions are met via the tableaus of the SRA method.\n\n\n\n"
},

{
    "location": "internals/solver_helpers.html#StochasticDiffEq.checkSRIOrder",
    "page": "Solver Extras",
    "title": "StochasticDiffEq.checkSRIOrder",
    "category": "Function",
    "text": "checkSRIOrder(RosslerSRI)\n\nDetermines whether the order conditions are met via the tableaus of the SRI method.\n\n\n\n"
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
    "location": "internals/tableaus.html#DiffEqDevTools.constructEuler",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructEuler",
    "category": "Function",
    "text": "Euler's method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRalston",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRalston",
    "category": "Function",
    "text": "Ralston's Order 2 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructHeun",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructHeun",
    "category": "Function",
    "text": "Heun's Order 2 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructKutta3",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructKutta3",
    "category": "Function",
    "text": "Kutta's Order 3 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructBS3",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructBS3",
    "category": "Function",
    "text": "constructBogakiShampine3()\n\nConstructs the tableau object for the Bogakai-Shampine Order 2/3 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructBogakiShampine3",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructBogakiShampine3",
    "category": "Function",
    "text": "constructBogakiShampine3()\n\nConstructs the tableau object for the Bogakai-Shampine Order 2/3 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRK4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRK4",
    "category": "Function",
    "text": "Classic RK4 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRK438Rule",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRK438Rule",
    "category": "Function",
    "text": "Classic RK4 3/8's rule method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRKF4",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRKF4",
    "category": "Function",
    "text": "Runge-Kutta-Fehberg Order 4/3\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRKF5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRKF5",
    "category": "Function",
    "text": "Runge-Kutta-Fehlberg Order 4/5 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructCashKarp",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructCashKarp",
    "category": "Function",
    "text": "constructCashKarp()\n\nConstructs the tableau object for the Cash-Karp Order 4/5 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructDormandPrince",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructDormandPrince",
    "category": "Function",
    "text": "constructDormandPrince()\n\nConstructs the tableau object for the Dormand-Prince Order 4/5 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructBS5",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructBS5",
    "category": "Function",
    "text": "An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine  Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructPapakostasPapaGeorgiou5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructPapakostasPapaGeorgiou5",
    "category": "Function",
    "text": "S.N. Papakostas and G. PapaGeorgiou higher error more stable\n\nA Family of Fifth-order Runge-Kutta Pairs, by S.N. Papakostas and G. PapaGeorgiou,  Mathematics of Computation,Volume 65, Number 215, July 1996, Pages 1165-1181.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructPapakostasPapaGeorgiou52",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructPapakostasPapaGeorgiou52",
    "category": "Function",
    "text": "S.N. Papakostas and G. PapaGeorgiou less stable lower error  Strictly better than DP5\n\nA Family of Fifth-order Runge-Kutta Pairs, by S.N. Papakostas and G. PapaGeorgiou,  Mathematics of Computation,Volume 65, Number 215, July 1996, Pages 1165-1181.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTsitouras5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTsitouras5",
    "category": "Function",
    "text": "Runge–Kutta pairs of orders 5(4) using the minimal set of simplifying assumptions,  by Ch. Tsitouras, TEI of Chalkis, Dept. of Applied Sciences, GR34400, Psahna, Greece.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLutherKonen5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLutherKonen5",
    "category": "Function",
    "text": "Luther and Konen's First Order 5 Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,  Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLutherKonen52",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLutherKonen52",
    "category": "Function",
    "text": "Luther and Konen's Second Order 5 Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,  Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLutherKonen53",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLutherKonen53",
    "category": "Function",
    "text": "Luther and Konen's Third Order 5 Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,  Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRungeFirst5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRungeFirst5",
    "category": "Function",
    "text": "Runge's First Order 5 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLawson5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLawson5",
    "category": "Function",
    "text": "Lawson's 5th order scheme\n\nAn Order Five Runge Kutta Process with Extended Region of Stability, J. Douglas Lawson,  Siam Journal on Numerical Analysis, Vol. 3, No. 4, (Dec., 1966) pages 593-597\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructSharpSmart5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructSharpSmart5",
    "category": "Function",
    "text": "Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum, by P.W.Sharp and E.Smart,  Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructBogakiShampine5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructBogakiShampine5",
    "category": "Function",
    "text": "An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine  Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructCassity5",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructCassity5",
    "category": "Function",
    "text": "Cassity's Order 5 method\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructButcher6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructButcher6",
    "category": "Function",
    "text": "Butcher's First Order 6 method\n\nOn Runge-Kutta Processes of High Order, by J. C. Butcher,  Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructButcher62",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructButcher62",
    "category": "Function",
    "text": "Butcher's Second Order 6 method\n\nOn Runge-Kutta Processes of High Order, by J. C. Butcher,  Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructButcher63",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructButcher63",
    "category": "Function",
    "text": "Butcher's Third Order 6\n\nOn Runge-Kutta Processes of High Order, by J. C. Butcher,  Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVernerRobust6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVernerRobust6",
    "category": "Function",
    "text": "From Verner's Website\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6A",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6A",
    "category": "Function",
    "text": "TanakaKasugaYamashitaYazaki Order 6 A\n\nOn the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,  by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,  Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6B",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6B",
    "category": "Function",
    "text": "constructTanakaKasugaYamashitaYazaki Order 6 B\n\nOn the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,  by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,  Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6C",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6C",
    "category": "Function",
    "text": "constructTanakaKasugaYamashitaYazaki Order 6 C\n\nOn the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,  by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,  Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6D",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTanakaKasugaYamashitaYazaki6D",
    "category": "Function",
    "text": "constructTanakaKasugaYamashitaYazaki Order 6 D\n\nOn the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,  by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,  Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructHuta6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructHuta6",
    "category": "Function",
    "text": "Anton Hutas First Order 6 method\n\nUne amélioration de la méthode de Runge-Kutta-Nyström pour la résolution numérique des équations différentielles du premièr ordre, by Anton Huta, Acta Fac. Nat. Univ. Comenian Math., Vol. 1, pages 201-224 (1956).\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructHuta62",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructHuta62",
    "category": "Function",
    "text": "Anton Hutas Second Order 6 method\n\nUne amélioration de la méthode de Runge-Kutta-Nyström pour la résolution numérique des équations différentielles du premièr ordre, by Anton Huta, Acta Fac. Nat. Univ. Comenian Math., Vol. 1, pages 201-224 (1956).\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVerner6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVerner6",
    "category": "Function",
    "text": "Verner Order 5/6 method\n\nA Contrast of a New RK56 pair with DP56, by Jim Verner,  Department of Mathematics. Simon Fraser University, Burnaby, Canada, 2006.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructDormandPrince6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructDormandPrince6",
    "category": "Function",
    "text": "Dormand-Prince Order 5//6 method\n\nP.J. Prince and J. R. Dormand, High order embedded Runge-Kutta formulae, Journal of Computational and Applied Mathematics . 7 (1981), pp. 67-75.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructSharpVerner6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructSharpVerner6",
    "category": "Function",
    "text": "Sharp-Verner Order 5/6 method\n\nCompletely Imbedded Runge-Kutta Pairs, by P. W. Sharp and J. H. Verner,  SIAM Journal on Numerical Analysis, Vol. 31, No. 4. (Aug., 1994), pages. 1169 to 1190.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructVern6",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructVern6",
    "category": "Function",
    "text": "From Verner's Website\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructClassicVerner6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructClassicVerner6",
    "category": "Function",
    "text": "EXPLICIT RUNGE-KUTFA METHODS WITH ESTIMATES OF THE LOCAL TRUNCATION ERROR\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructChummund6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructChummund6",
    "category": "Function",
    "text": "Chummund's First Order 6 method\n\nA three-dimensional family of seven-step Runge-Kutta methods of order 6, by G. M. Chammud (Hammud), Numerical Methods and programming, 2001, Vol.2, 2001, pages 159-166 (Advanced Computing Scientific journal published by the Research Computing Center of the Lomonosov Moscow State Univeristy)\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructChummund62",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructChummund62",
    "category": "Function",
    "text": "Chummund's Second Order 6 method\n\nA three-dimensional family of seven-step Runge-Kutta methods of order 6, by G. M. Chammud (Hammud), Numerical Methods and programming, 2001, Vol.2, 2001, pages 159-166 (Advanced Computing Scientific journal published by the Research Computing Center of the Lomonosov Moscow State Univeristy)\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructPapakostas6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructPapakostas6",
    "category": "Function",
    "text": "Papakostas's Order 6\n\nOn Phase-Fitted modified Runge-Kutta Pairs of order 6(5), by Ch. Tsitouras and I. Th. Famelis,  International Conference of Numerical Analysis and Applied Mathematics, Crete, (2006)\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructLawson6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructLawson6",
    "category": "Function",
    "text": "Lawson's Order 6\n\nAn Order 6 Runge-Kutta Process with an Extended Region of Stability, by J. D. Lawson,  Siam Journal on Numerical Analysis, Vol. 4, No. 4 (Dec. 1967) pages 620-625.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTsitourasPapakostas6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTsitourasPapakostas6",
    "category": "Function",
    "text": "Tsitouras-Papakostas's Order 6\n\nCheap Error Estimation for Runge-Kutta methods, by Ch. Tsitouras and S.N. Papakostas, Siam Journal on Scientific Computing, Vol. 20, Issue 6, Nov 1999.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructDormandLockyerMcCorriganPrince6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructDormandLockyerMcCorriganPrince6",
    "category": "Function",
    "text": "DormandLockyerMcCorriganPrince Order 6 Global Error Estimation\n\nGlobal Error estimation with Runge-Kutta triples, by J.R.Dormand, M.A.Lockyer, N.E.McCorrigan and P.J.Prince,  Computers and Mathematics with Applications, 18 (1989) pages 835-846.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVernerEfficient6",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVernerEfficient6",
    "category": "Function",
    "text": "From Verner's Website\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructMikkawyEisa",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructMikkawyEisa",
    "category": "Function",
    "text": "Mikkawy-Eisa Order 6\n\nA general four-parameter non-FSAL embedded Runge–Kutta algorithm of orders 6 and 4 in seven stages,  by M.E.A. El-Mikkawy and M.M.M. Eisa,  Applied Mathematics and Computation, Vol. 143, No. 2, (2003) pages 259 to 267.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVernerEfficient7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVernerEfficient7",
    "category": "Function",
    "text": "From Verner's website\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructClassicVerner7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructClassicVerner7",
    "category": "Function",
    "text": "EXPLICIT RUNGE-KUTFA METHODS WITH ESTIMATES OF THE LOCAL TRUNCATION ERROR\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructSharpVerner7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructSharpVerner7",
    "category": "Function",
    "text": "Completely Imbedded Runge-Kutta Pairs, by P.W.Sharp and J.H.Verner, Siam Journal on Numerical Analysis, Vol.31, No.4. (August 1994) pages 1169-1190.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTanakaYamashitaStable7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTanakaYamashitaStable7",
    "category": "Function",
    "text": "On the Optimization of Some Nine-Stage Seventh-order Runge-Kutta Method, by M. Tanaka, S. Muramatsu and S. Yamashita, Information Processing Society of Japan, Vol. 33, No. 12 (1992) pages 1512-1526.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructSharpSmart7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructSharpSmart7",
    "category": "Function",
    "text": "Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum, by P.W.Sharp and E.Smart,  Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTanakaYamashitaEfficient7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTanakaYamashitaEfficient7",
    "category": "Function",
    "text": "On the Optimization of Some Nine-Stage Seventh-order Runge-Kutta Method, by M. Tanaka, S. Muramatsu and S. Yamashita, Information Processing Society of Japan, Vol. 33, No. 12 (1992) pages 1512-1526.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVernerRobust7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVernerRobust7",
    "category": "Function",
    "text": "From Verner's website\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructTanYam7",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructTanYam7",
    "category": "Function",
    "text": "On the Optimization of Some Nine-Stage Seventh-order Runge-Kutta Method, by M. Tanaka, S. Muramatsu and S. Yamashita, Information Processing Society of Japan, Vol. 33, No. 12 (1992) pages 1512-1526.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructEnrightVerner7",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructEnrightVerner7",
    "category": "Function",
    "text": "The Relative Efficiency of Alternative Defect Control Schemes for High-Order Continuous Runge-Kutta Formulas  W. H. Enright SIAM Journal on Numerical Analysis, Vol. 30, No. 5. (Oct., 1993), pp. 1419-1445.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructDormandPrince8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructDormandPrince8",
    "category": "Function",
    "text": "constructDormandPrice8()\n\nConstructs the tableau object for the Dormand-Prince Order 6/8 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructRKF8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructRKF8",
    "category": "Function",
    "text": "constructRKF8()\n\nConstructs the tableau object for the Runge-Kutta-Fehlberg Order 7/8 method.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructCooperVerner8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructCooperVerner8",
    "category": "Function",
    "text": "Some Explicit Runge-Kutta Methods of High Order, by G. J. Cooper and J. H. Verner,  SIAM Journal on Numerical Analysis, Vol. 9, No. 3, (September 1972), pages 389 to 405\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructCooperVerner82",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructCooperVerner82",
    "category": "Function",
    "text": "Some Explicit Runge-Kutta Methods of High Order, by G. J. Cooper and J. H. Verner,  SIAM Journal on Numerical Analysis, Vol. 9, No. 3, (September 1972), pages 389 to 405\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTsitourasPapakostas8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTsitourasPapakostas8",
    "category": "Function",
    "text": "Cheap Error Estimation for Runge-Kutta methods, by Ch. Tsitouras and S.N. Papakostas,  Siam Journal on Scientific Computing, Vol. 20, Issue 6, Nov 1999.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructEnrightVerner8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructEnrightVerner8",
    "category": "Function",
    "text": "The Relative Efficiency of Alternative Defect Control Schemes for High-Order Continuous Runge-Kutta Formulas  W. H. Enright SIAM Journal on Numerical Analysis, Vol. 30, No. 5. (Oct., 1993), pp. 1419-1445.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructdverk78",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructdverk78",
    "category": "Function",
    "text": "Jim Verner's \"Maple\" (dverk78)\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructClassicVerner8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructClassicVerner8",
    "category": "Function",
    "text": "EXPLICIT RUNGE-KUTFA METHODS WITH ESTIMATES OF THE LOCAL TRUNCATION ERROR\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructDormandPrince8_64bit",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructDormandPrince8_64bit",
    "category": "Function",
    "text": "constructDormandPrice8_64bit()\n\nConstructs the tableau object for the Dormand-Prince Order 6/8 method with the approximated coefficients from the paper. This works until below 64-bit precision.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructCurtis8",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructCurtis8",
    "category": "Function",
    "text": "An Eighth Order Runge-Kutta process with Eleven Function Evaluations per Step, by A. R. Curtis,  Numerische Mathematik, Vol. 16, No. 3 (1970), pages 268 to 277\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructTsitPap8",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructTsitPap8",
    "category": "Function",
    "text": "Cheap Error Estimation for Runge-Kutta methods, by Ch. Tsitouras and S.N. Papakostas,  Siam Journal on Scientific Computing, Vol. 20, Issue 6, Nov 1999.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructSharp9",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructSharp9",
    "category": "Function",
    "text": "Journal of Applied Mathematics & Decision Sciences, 4(2), 183-192 (2000),  \"High order explicit Runge-Kutta pairs for ephemerides of the Solar System and the Moon\".\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTsitouras9",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTsitouras9",
    "category": "Function",
    "text": "Optimized explicit Runge-Kutta pairs of order 9(8), by Ch. Tsitouras,  Applied Numerical Mathematics, 38 (2001) 123-134.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructTsitouras92",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructTsitouras92",
    "category": "Function",
    "text": "Optimized explicit Runge-Kutta pairs of order 9(8), by Ch. Tsitouras,  Applied Numerical Mathematics, 38 (2001) 123-134.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVernerEfficient9",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVernerEfficient9",
    "category": "Function",
    "text": "From Verner's Webiste\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructVern9",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructVern9",
    "category": "Function",
    "text": "From Verner's Webiste\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVerner916",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVerner916",
    "category": "Function",
    "text": "Verner 1991 First Order 5/6 method\n\nSome Ruge-Kutta Formula Pairs, by J.H.Verner,  SIAM Journal on Numerical Analysis, Vol. 28, No. 2 (April 1991), pages 496 to 511.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVerner9162",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVerner9162",
    "category": "Function",
    "text": "Verner 1991 Second Order 5/6 method\n\nSome Ruge-Kutta Formula Pairs, by J.H.Verner,  SIAM Journal on Numerical Analysis, Vol. 28, No. 2 (April 1991), pages 496 to 511.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructVernerRobust9",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructVernerRobust9",
    "category": "Function",
    "text": "From Verner's Webiste\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructFeagin10",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructFeagin10",
    "category": "Function",
    "text": "constructFeagin10\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructFeagin10Tableau",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructFeagin10Tableau",
    "category": "Function",
    "text": "Feagin10 in Tableau form\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructOno10",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructOno10",
    "category": "Function",
    "text": "Ono10\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructCurtis10",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructCurtis10",
    "category": "Function",
    "text": "High-order Explicit Runge-Kutta Formulae, Their uses, and Limitations, A.R.Curtis, J. Inst. Maths Applics (1975) 16, 35-55.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructHairer10",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructHairer10",
    "category": "Function",
    "text": "A Runge-Kutta Method of Order 10, E. Hairer, J. Inst. Maths Applics (1978) 21, 47-59.\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructBaker10",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructBaker10",
    "category": "Function",
    "text": "Tom Baker, University of Teeside. Part of RK-Aid http://www.scm.tees.ac.uk/users/u0000251/research/researcht.htm http://www.scm.tees.ac.uk/users/u0000251/j.r.dormand/t.baker/rk10921m/rk10921m\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructFeagin12",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructFeagin12",
    "category": "Function",
    "text": "constructFeagin12\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructOno12",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructOno12",
    "category": "Function",
    "text": "On the 25 stage 12th order explicit Runge-Kutta method, by Hiroshi Ono. Transactions of the Japan Society for Industrial and applied Mathematics, Vol. 6, No. 3, (2006) pages 177 to 186\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructFeagin12Tableau",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructFeagin12Tableau",
    "category": "Function",
    "text": "Tableau form of Feagin12\n\n\n\n"
},

{
    "location": "internals/tableaus.html#OrdinaryDiffEq.constructFeagin14",
    "page": "ODE Tableaus",
    "title": "OrdinaryDiffEq.constructFeagin14",
    "category": "Function",
    "text": "constructFeagin14\n\n\n\n"
},

{
    "location": "internals/tableaus.html#DiffEqDevTools.constructFeagin14Tableau",
    "page": "ODE Tableaus",
    "title": "DiffEqDevTools.constructFeagin14Tableau",
    "category": "Function",
    "text": "Tableau form of Feagin14\n\n\n\n"
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
