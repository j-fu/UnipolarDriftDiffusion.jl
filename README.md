# Introduction

`UnipolarDriftDiffusion.jl` provides numerical examples for finite volume schemes
 for the paper "A numerical analysis focused comparison of several Finite Volume schemes for an
 Unipolar Degenerated Drift-Diffusion Model"
 by C.  Cancès, C. Chainais-Hilliaret, J. Fuhrmann, B. Gaudeul (submitted)

It uses the finite volume solver [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)

## Usage

For generating the  figures, perform
````julia
repl> using Pkg
repl> Pkg.activate(".")
repl> ENV["MPLBACKEND"]="agg"
repl> include("UnipolarDriftDiffusion.jl")
repl> UnipolarDriftDiffusion.makefigures(quick=true)
````

`quick=false` generates figures as they are in the paper, but results in an overnight calculation.

# Module header and packges used

```julia
module UnipolarDriftDiffusion
```

## Standard Julia modules

```julia
using Printf
using PyPlot
using LinearAlgebra
```

## Generate documentation output from Julia code

```julia
using Literate
```

## Binary file input/output

```julia
using FileIO
```

## Finite volume solver

```julia
using VoronoiFVM
```

## Resizable multidimensional arrays

```julia
using ElasticArrays
```

# Flux and reaction callbacks and their data


## Data provided to callback functions

```julia
mutable struct Data <: VoronoiFVM.AbstractData
    iΦ::Int32 # Potential index in species list
    ic::Int32  # Concentration index in species list
    doping::Float64 # Doping
    D::Float64      # Mobility coefficient
    Data()=new()
end
```

## Sedan flux

Appearantly first describe by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

see also the 198? Fortran code available via http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

```julia
function sflux!(f,u,edge,data)

    uk=viewK(edge,u) # unknowns at left end of discretization edge
    ul=viewL(edge,u) # unknowns at right end of discretization edge
    iΦ=data.iΦ # index of potential
    ic=data.ic  # index of concentration

    # Poisson flux
    f[iΦ]=uk[iΦ]-ul[iΦ]

    # Caclculate excess chemical potentials
    muk=-log(1-uk[ic])
    mul=-log(1-ul[ic])

    # Combine potential gradient and excess chemical gradient
    arg=(uk[iΦ]-ul[iΦ])+(muk-mul)

    # Call Bernoulli function
    bp,bm=fbernoulli_pm(arg)

    # Calculate drift-diffusion flux
    f[ic]=data.D*(bm*uk[ic]-bp*ul[ic])

end
```

## Flux based on averaging activity coefficients

J. Fuhrmann, DOI 10.1016/j.cpc.2015.06.004

```julia
function aflux!(f,u,edge,data)
    uk=viewK(edge,u)
    ul=viewL(edge,u)
    iΦ=data.iΦ
    ic=data.ic

    f[iΦ]=uk[iΦ]-ul[iΦ]
    ck=uk[ic]
    cl=ul[ic]

    # Calculate activities
    ak=ck/(1.0-ck)
    al=cl/(1.0-cl)

    # Keep only potential difference in Bernoulli argument
    arg=uk[iΦ]-ul[iΦ]
    bp,bm=fbernoulli_pm(arg)
```

Multiply weighted finite difference of  activties with
averaged inverse activity coefficients

```julia
    f[ic]=data.D*(0.5*(1.0/(1.0+ak)+ 1.0/(1.0+al))*(bm*ak-bp*al))
end
```

## Central difference flux

B. Andreianov, C. Cancès, A. Moussa, DOI  10.1016/j.jfa.2017.08.010

```julia
function cflux!(f,u,edge,data)
    uk=viewK(edge,u)
    ul=viewL(edge,u)
    iΦ=data.iΦ
    ic=data.ic
    f[iΦ]=(uk[iΦ]-ul[iΦ])
    ck=uk[ic]
    cl=ul[ic]
    hk=log(ck/(1.0-ck))
    hl=log(cl/(1.0-cl))
    f[ic]=data.D*(0.5*(ck+cl)*(hk-hl+(uk[iΦ]-ul[iΦ])))
end
```

## Bessemoulin-Chatard flux

M. Bessemoulin-Chatard, DOI 10.1007/s00211-012-0448-x

```julia
function bflux!(f,u,edge,data)
    uk=viewK(edge,u)
    ul=viewL(edge,u)
    iΦ=data.iΦ
    ic=data.ic
    f[iΦ]=(uk[iΦ]-ul[iΦ])
    ck=uk[ic]
    cl=ul[ic]
    gk=1.0/(1.0-ck)
    gl=1.0/(1.0-cl)

    if isapprox(ck,cl,atol=1.0e-10)
        dkl=0.5*(gk+gl)
    else
        dkl=(log(gk*ck)-log(gl*cl))/(log(ck)-log(cl))
    end
    arg=(uk[iΦ]-ul[iΦ])/dkl
    bp,bm=fbernoulli_pm(arg)
    f[ic]=data.D*(dkl*(bm*ck-bp*cl))
end
```

## Storage term

 This describes the term under the time derivative

```julia
function storage!(f,u,node,data)
    # Poisson equation is stationary
    f[data.iΦ]=0
    # Concentration
    f[data.ic]=u[data.ic]
end
```

## Reaction term

 This describes charge density and  recombination

```julia
function reaction!(f,u,node,data)
    # Charge density
    f[data.iΦ]=-data.doping-u[data.ic]
    # Zero recombination/generation
    f[data.ic]=0
end
```

## Combine physical data of the problem

```julia
function DefaultPhysics(;flux=sflux!)
    # Data record for function callbacs
    data=Data()
    data.iΦ=1
    data.ic=2
    data.doping=-0.5
    data.D=1

    # Physics combines data and callback function
    physics=VoronoiFVM.Physics(data=data,
                               num_species=2,
                               flux=flux,
                               reaction=reaction!,
                               storage=storage!
                               )
    return physics
end
```

# Simulation examples


## Evolution of constant concentrations into equilibrium with different applied potentials

```julia
function run_to_equilibrium(;
                            n=100,        # number of nodes
                            tend=1.0e4,   # end of time evolution
                            pyplot=false, # create plots during simulation
                            genplot=false,# generate plots for paper (*eevol.pdf, *tevol.pdf)
                            subcase="0",  # subcases ("1","2","3" for paper)
                            cases=["s"])  # Fluxes to consider

    # Create discretization grid
    L=50
    h=L/convert(Float64,n)
    X=collect(0:h:L)
    grid=VoronoiFVM.Grid(X)

    # Create dict where to put results
    alldata=DataDict()

    # Declare some  variables used later
    ic=nothing
    iΦ=nothing
    data=nothing
    U=nothing
    sys=nothing

    # Loop over different fluxes
    for case in cases
        physics=DefaultPhysics(flux=alldata[case]["flux"])
        println(physics.flux)
        data=physics.data
        iΦ=data.iΦ
        ic=data.ic

        # Create finite volume system
        sys=VoronoiFVM.DenseSystem(grid,physics)
        enable_species!(sys,1,[1])
        enable_species!(sys,2,[1])

        # Declare Dirichlet bc
        sys.boundary_factors[data.iΦ,1]=VoronoiFVM.Dirichlet
        sys.boundary_factors[data.iΦ,2]=VoronoiFVM.Dirichlet


        # Set up testfunction fot flux calculations
        factory=VoronoiFVM.TestFunctionFactory(sys)
        tf=testfunction(factory,[2],[1])

        # Create initial values
        inival=unknowns(sys)

        # Default case
        @views inival[data.iΦ,:].=0
        @views inival[data.ic,:].=0.7
        Φ_0=0
        sys.boundary_values[data.iΦ,1]=Φ_0
        sys.boundary_values[data.iΦ,2]=0.0


        # Potential jump
        if subcase=="1"
            @views inival[data.ic,:].=0.5
            Φ_0=10
            sys.boundary_values[data.iΦ,1]=Φ_0
            sys.boundary_values[data.iΦ,2]=0.0
        end

        # Depleted cell
        if subcase=="2"
            @views inival[data.ic,:].=0.3
            Φ_0=0
            sys.boundary_values[data.iΦ,1]=Φ_0
            sys.boundary_values[data.iΦ,2]=0.0
        end

        # Oversaturated cell
        if subcase=="3"
            @views inival[data.ic,:].=0.7
            Φ_0=0
            sys.boundary_values[data.iΦ,1]=Φ_0
            sys.boundary_values[data.iΦ,2]=0.0
        end



        # Stationary solution of Poisson equation
        # Immobilise c, one large timestep
        U=unknowns(sys)
        data.D=0
        solve!(U,inival,sys,tstep=1.0e30)
        data.D=1

        inival.=U
        if pyplot
            plot_solution(sys,U)
        end


        # Control data for damped newton
        control=VoronoiFVM.NewtonControl()
        control.verbose=false
        control.max_lureuse=0
        control.damp_initial=0.01

        # Set up result arrays
        times=zeros(0)
        Energies=zeros(0)
        Φ=ElasticArray{Float64}(undef,length(X),1)
        C=ElasticArray{Float64}(undef,length(X),1)

        # Push initial value
        Φ[:,1].=inival[iΦ,:]
        C[:,1].=inival[ic,:]
        push!(alldata[case],Pair("Φ",Φ))
        push!(alldata[case],Pair("C",C))
        push!(alldata[case],Pair("times",times))
        push!(alldata[case],Pair("Energies",Energies))

        # Time evolution
        t=0.0
        tstep=1.0e-4
        dtgrowth=1.15
        istep=0
        while t<tend
            # Solve for new timestep
            t+=tstep
            istep+=1
            solve!(U,inival,sys,control=control,tstep=tstep)

            # Calculate current (aka \nabla \phi_0)
            I=integrate(sys,tf,U)
            gradΦ_0=I[iΦ]

            # Calculate free energy
            E_Φ=-Φ_0*gradΦ_0
            E_c=0
            H= (c) -> c*log(c) + (1-c)*log(1-c)
            for i=1:length(X)-1
                h=(X[i+1]-X[i])
                gradΦ=(U[iΦ,i+1]-U[iΦ,i])/h
                E_Φ+=0.5*h*gradΦ^2
                E_c+=0.5*h*(H(U[ic,i+1])+H(U[ic,i]))
            end

            @printf("time=%8.4e gradΦ_0=%e E=%e\n",t, gradΦ_0,E_Φ+E_c)

            # Push results
            push!(times,t)
            push!(Energies,E_Φ+E_c)
            resize!(Φ,length(X),istep)
            resize!(C,length(X),istep)
            Φ[:,istep].=U[iΦ,:]
            C[:,istep].=U[ic,:]

            if pyplot
                plot_solution(sys,U)
            end

            # Prepare next timestep
            tstep*=dtgrowth
            control.damp_initial=1
            inival.=U
        end
    end


    if genplot

        # Log-log plot of relative free energy
        figure(num=1)
        PyPlot.clf()
        for case in cases
            times=alldata[case]["times"]
            Energies=alldata[case]["Energies"]
            line=alldata[case]["line"]
            flux=alldata[case]["name"]
            marker=alldata[case]["marker"]
            markersize=alldata[case]["markersize"]
            PyPlot.loglog(times,Energies.+(1.0e-13-Energies[end]),line,label=flux, marker=marker,markersize=markersize)
        end
        PyPlot.xlabel("\$t\$")
        PyPlot.grid()
        PyPlot.legend(loc="lower left")
        PyPlot.pause(1.0e-10)

        fname="eevol-$(subcase).pdf"
        mysavefig(fname)
        PyPlot.close() # need to close here as clf() seems to keep some of the state...


        # Semilog-t plot of time evolution
        figure(num=1,figsize=[6,4])
        PyPlot.clf()
        case="s"
        times=alldata[case]["times"]
        Φ=alldata[case]["Φ"]
        C=alldata[case]["C"]
        Φmax=max(maximum(Φ),-minimum(Φ))
        allsol=Array{Float64}(undef,2,length(X)*length(times))
        allsol[ic,:]=C[:]
        allsol[iΦ,:]=Φ[:]
        plot2d(data,X,times,allsol,umax=Φmax,uraster=5,ysc="log",ylbl="t")
        PyPlot.tight_layout()
        fname="tevol-$(subcase).pdf"
        mysavefig(fname)
        PyPlot.close() # need to close here as clf() seems to keep some of the state...
    end
return U,sys
end
```

## Stationary solution

```julia
function run_stationary(;L=50,       # Domain size
                        case="s",    # Case to consider
                        n=100,       # Number  of unknowns
                        pyplot=false # Plot during calculation
                        )
    # Initialize data, grid
    alldata=DataDict()
    physics=DefaultPhysics(flux=alldata[case]["flux"])
    iΦ=physics.data.iΦ
    ic=physics.data.ic
    h=L/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:L))

    # Create finite volume system
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])

    # Set boundary conditions
    sys.boundary_values[iΦ,1]=0
    sys.boundary_values[iΦ,2]=0

    sys.boundary_factors[iΦ,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[iΦ,2]=VoronoiFVM.Dirichlet

    sys.boundary_values[ic,1]=0.5
    sys.boundary_factors[ic,1]=VoronoiFVM.Dirichlet

    sys.boundary_values[ic,2]=0.5
    sys.boundary_factors[ic,2]=VoronoiFVM.Dirichlet

    # Initial values
    inival=unknowns(sys)
    @views inival[iΦ,:].=0.5
    @views inival[ic,:].=0.5

    control=VoronoiFVM.NewtonControl()
    control.verbose=true
    control.damp_initial=0.01
    control.damp_growth=2
    control.tol_absolute=1.0e-12
    control.Δp=0.1
    control.Δp_min=1.0e-5
    control.Δp_max=0.1

    delta=1.0e-4
    sol=unknowns(sys)
    # Pre-step function for embedding (to critical boundary concentraions)
    pre=function(sol,p)
        px=p^0.25
        bc1=0.9+px*(0.999-0.9)
        bc2=1.0-bc1
        @printf("n=%d px=%e bc1=%e bc2=%e\n",n,px,bc1,bc2)
        sys.boundary_values[ic,2]=bc1
        sys.boundary_values[ic,1]=bc2
    end

    # Post-step function for embedding
    post=function(sol,p)
        if pyplot
            plot_solution(sys,sol)
        end
    end
    # Obtain stationary solution via embedding
    @time embed!(sol,inival,sys, control=control,pre=pre, post=post)
    return sol,sys
end
```

## Compare one-dimensional solutions (either stationary, or transient moment of time)

Generates *-schemes*.pdf, *-refsol.pdf

```julia
function run_compare(;maxref=10,stationary=true)
    iΦ=1
    ic=2
    nmax=10*2^maxref
    L=50
    @printf("nmax=%d\n",nmax)
    @printf("hmin=%e\n",1.0/convert(Float64,nmax))

    # Create reference solution
    if stationary
        pfx="stat-"
        refsol,refsys=run_stationary(L=L,case="s",n=nmax,pyplot=false)
    else
        pfx="tran-"
        refsol,refsys=run_to_equilibrium(n=nmax,cases=["s"],tend=1.0e1, subcase="1")
    end

    write("refsol.dat",refsol)
    plot_solution(refsys,refsol)
    PyPlot.tight_layout()
    mysavefig("$(pfx)refsol.pdf")


    hmin=L/convert(Float64,nmax)
    alldata=DataDict()

    # Loop over fluxes
    cases=["c","a","b","s"]
    for case in cases
        xflux=alldata[case]["flux"]
        flux=alldata[case]["name"]
        l2line=alldata[case]["line"]

        # Calculate gradient of reference solution
        refflux=zeros(2,nmax)
        for i=1:nmax
            refflux[iΦ,i]=(refsol[iΦ,i+1]-refsol[iΦ,i])/hmin
            refflux[ic,i]=(refsol[ic,i+1]-refsol[ic,i])/hmin
        end


        H=zeros(0)
        L2Error=zeros(0)
        H1Error=zeros(0)

        push!(alldata[case],Pair("H",H))
        push!(alldata[case],Pair("L2Error",L2Error))
        push!(alldata[case],Pair("H1Error",H1Error))

        # Loop over refinement levels
        xmaxref=maxref-2
        for ref=1:xmaxref
            nn=10*2^ref
            h=1.0/convert(Float64,nn)
            # Obtain stationary resp. transient solution
            if stationary
                sol,sys=run_stationary(L=L,case=case,n=nn)
            else
                sol,sys=run_to_equilibrium(n=nn,cases=[case],tend=1.0e1, subcase="1")
            end

            # Interpolate solution to finest level of reference solution
            intersol=copy(refsol)
            intersol.=0
            ninter=2^(maxref-ref)
            imax=1
            for i=1:nn
                phi1=sol[iΦ,i]
                phi2=sol[iΦ,i+1]
                c1=sol[ic,i]
                c2=sol[ic,i+1]
                for iinter=1:ninter
                    x2=(iinter-1)/ninter
                    x1=1.0-x2
                    intersol[iΦ,imax]=phi1*x1+phi2*x2
                    intersol[ic,imax]=c1*x1+c2*x2
                    imax=imax+1
                end
            end
            intersol[iΦ,nmax+1]=sol[iΦ,nn+1]
            intersol[ic,nmax+1]=sol[ic,nn+1]

            # Calculate gradient of interpolated solution
            interflux=zeros(2,nmax)
            for i=1:nmax
                interflux[iΦ,i]=(intersol[iΦ,i+1]-intersol[iΦ,i])/hmin
                interflux[ic,i]=(intersol[ic,i+1]-intersol[ic,i])/hmin
            end

            # Calculate norms
            l2error=sqrt(hmin)*norm(refsol-intersol,2)
            h1error=sqrt(hmin)*norm(refflux-interflux,2)
            @printf("l2 error: %e\n",l2error)
            @printf("h1 error: %e\n",h1error)

            # Push result
            push!(H,h)
            push!(L2Error,l2error)
            push!(H1Error,h1error)
        end
        print(H)
        print(L2Error)
        print(H1Error)
    end

    # Generate pdfs for plots.
    PyPlot.clf()
    H=nothing
    for case in cases
        L2Error=alldata[case]["L2Error"]
        H=alldata[case]["H"]
        line=alldata[case]["line"]
        flux=alldata[case]["name"]
        PyPlot.loglog(H,L2Error,line, label=flux, linewidth=2.5, marker="o",markersize=5)
    end
    PyPlot.loglog(H,4.0e2*H.*H,"k--",label="\$O(h^2)\$")
    PyPlot.xlabel("\$h\$")
    PyPlot.ylabel("\$L^2\$ Error")
    PyPlot.legend(loc="lower right")
    PyPlot.grid()
    PyPlot.xlim(5.0e-5,1.0e-1)
    PyPlot.tight_layout()
    mysavefig("$(pfx)schemes-l2.pdf")


    PyPlot.clf()
    H=nothing
    for case in cases
        H1Error=alldata[case]["H1Error"]
        H=alldata[case]["H"]
        line=alldata[case]["line"]
        flux=alldata[case]["name"]
        marker=alldata[case]["marker"]
        markersize=alldata[case]["markersize"]
        PyPlot.loglog(H,H1Error,line, label=flux, linewidth=2.5, marker=marker,markersize=markersize)
    end
    PyPlot.loglog(H,0.5e2*H,"k--",label="\$O(h)\$")
    PyPlot.xlabel("\$h\$")
    PyPlot.ylabel("\$H^1\$ Error")
    PyPlot.legend(loc="lower right")
    PyPlot.grid()
    PyPlot.xlim(5.0e-5,1.0e-1)
    PyPlot.tight_layout()
    mysavefig("$(pfx)schemes-h1.pdf")
    PyPlot.close()

end
```

## Bias loop for field effect transistor

```julia
function run_fet(;nref=0,
                 plotgrid=false,
                 postplot=false,
                 genref=false,
                 genplot=false,
                 clearplot=true,
                 conv=false)

    L=0.1
    wcontact=0.2*L
    wgate=0.4*L
    wy=0.1*L
    nx=10*2^nref
    ny=5*2^nref

    ugmax=50
    udmax=20

    uplotmax=25

    igate=1
    isource=2
    idrain=3
    ineutral=4
    d=0.1*wy
    pen=1.0e10

    cases=["c","a","b","s"]
    X=collect(0.0:L/nx:L)
    Y=collect(0.0:wy/ny:wy)

    grid=VoronoiFVM.Grid(X,Y)
    bfacemask!(grid,[0,0],[L,wy],ineutral)
    bfacemask!(grid,[0,wy],[wcontact,wy],isource)
    bfacemask!(grid,[L-wcontact,wy],[L,wy],idrain)
    bfacemask!(grid,[0.5*(L-wgate),wy],[L/2+wgate/2,wy],igate)

    doping=zeros(num_nodes(grid))

    @printf("num_nodes=%d\n", num_nodes(grid))


    if plotgrid
        figure(num=1,figsize=[5.0,1.5])
        clf()
        fvmplot(grid)
        show()
        return
    end


    ugate=0.0


    breaction=function(f,u,node,data)
        if  node.region==igate
            f[data.iΦ]=(u[data.iΦ]-ugate)/d
            f[data.ic]=0
        else
            f[1]=0
            f[2]=0
        end
    end


    data=Data()
    data.iΦ=1
    data.ic=2
    data.doping=-0.5
    data.D=1.0

    alldata=DataDict()

    if genref
        cases=["s"]
    end
    if genplot
        cases=["s"]
    end

    for case in cases
        xflux=alldata[case]["flux"]
        flux=alldata[case]["name"]
        l2line=alldata[case]["line"]


        physics=VoronoiFVM.Physics(data=data,
                                   num_species=2,
                                   flux= xflux,
                                   reaction=reaction!,
                                   storage=storage!,
                                   breaction=breaction
                                   )

        println("sys:")
        @time sys=VoronoiFVM.DenseSystem(grid,physics)
        println("enable spec:")
        @time enable_species!(sys,1,[1])
        println("enable spec:")
        @time enable_species!(sys,2,[1])


        println("factory:")
        @time factory=VoronoiFVM.TestFunctionFactory(sys)
        println("tf:")
        @time tf_drain=testfunction(factory,[isource],[idrain])



        sys.boundary_values[data.ic,idrain]=0.5
        sys.boundary_factors[data.ic,idrain]=VoronoiFVM.Dirichlet

        sys.boundary_values[data.ic,isource]=0.5
        sys.boundary_factors[data.ic,isource]=VoronoiFVM.Dirichlet

        sys.boundary_values[data.iΦ,idrain]=0
        sys.boundary_factors[data.iΦ,idrain]=VoronoiFVM.Dirichlet

        sys.boundary_values[data.iΦ,isource]=0
        sys.boundary_factors[data.iΦ,isource]=VoronoiFVM.Dirichlet

        sys.boundary_values[data.iΦ,igate]=0
        sys.boundary_factors[data.iΦ,igate]=1.0/d


        inival=unknowns(sys)
        @views inival[data.iΦ,:].=0.0
        @views inival[data.ic,:].=0.5


        sol=unknowns(sys)
        ugate=0.0

        control=VoronoiFVM.NewtonControl()
        control.damp_initial=0.5
        control.damp_growth=2
        control.Δp=0.2
        control.Δp_min=1.0e-5
        control.Δp_max=0.5
        control.verbose=true
        control.tol_relative=1.0e-10
        control.tol_absolute=1.0e-10
        control.tol_round=1.0e-10
        control.tol_linear=1.0e-4
        control.max_round=4
        control.max_iterations=20
        control.max_lureuse=0
        control.handle_exceptions=true

        print(control)

        @time solve!(sol,inival,sys, control=control)


        inival.=sol
        function pre(sol,p)
            ugate=p*ugmax
            sys.boundary_values[data.iΦ,isource]=p*udmax/2
            sys.boundary_values[data.iΦ,idrain]=-p*udmax/2
        end

        @time embed!(sol,inival,sys, control=control,pre=pre)

        inival.=sol


        Idrain=zeros(0)
        Ugate=collect(ugmax:-0.1*ugmax:-ugmax)

        push!(alldata[case],Pair("Ugate",Ugate))
        push!(alldata[case],Pair("Idrain",Idrain))

        I=integrate(sys,tf_drain,sol)
        @printf("u_gate=%e I_drain=%e\n",ugate,I[data.ic])
        push!(Idrain,I[data.ic])
        if postplot
            plot2d(data,X,Y,sol,umax=uplotmax)
        end
        if genplot
            figure(num=1,figsize=[15,3])
            plot2d(data,X,Y,sol,umax=uplotmax)
            mysavefig("fet-closed.pdf")
            close()
        end


        for i=1:length(Ugate)-1
            embed!(sol,inival,sys,control=control,
                   pre=function(sol,p)
                   ugate=Ugate[i]+p*(Ugate[i+1]-Ugate[i])
                   end,
                   )
            I=integrate(sys,tf_drain,sol)
            @printf("u_gate=%e  I_drain=%e\n",ugate,I[data.ic])
            @assert(ugate≈Ugate[i+1])
            push!(Idrain,I[data.ic])
            if postplot
                plot2d(data,X,Y,sol,umax=uplotmax)
            end
            if genplot
                if isapprox(ugate,0.0,atol=1.0e-2)
                    figure(num=1,figsize=[15,3])
                    plot2d(data,X,Y,sol,umax=uplotmax)
                    mysavefig("fet-g0.pdf")
                    close()
                end
            end

            inival.=sol
        end
        if genref
            write("fet-refsol.dat",Idrain)
        end
        if genplot
            figure(num=1,figsize=[15,3])
            plot2d(data,X,Y,sol,umax=uplotmax)
            mysavefig("fet-opened.pdf")
            close()
        end
    end
    if conv
       return alldata
    end

if genplot
    PyPlot.clf()
    for case in cases
        Idrain=alldata[case]["Idrain"]
        Ugate=alldata[case]["Ugate"]
        line=alldata[case]["line"]
        flux=alldata[case]["name"]
        marker=alldata[case]["marker"]
        markersize=alldata[case]["markersize"]
        PyPlot.plot(Ugate,-Idrain,line,label=flux, marker=marker,markersize=markersize)
    end
    PyPlot.grid()
    PyPlot.legend(loc="lower left")

    show()
    mysavefig(@sprintf("fet-iv-nref%d.pdf",nref))
    close()
end


end
```

## Convergence test for unipolar field effect transistor

```julia
function run_conv(;maxref=3,genref=false, genconv=true)

    cases=["c","a","b","s"]
    if genconv
        refcurr=zeros(21)
        read!("fet-refsol.dat",refcurr)

        allres=nothing
        alldata=DataDict()

        push!(alldata["a"],Pair("error",zeros(0)))
        push!(alldata["b"],Pair("error",zeros(0)))
        push!(alldata["c"],Pair("error",zeros(0)))
        push!(alldata["s"],Pair("error",zeros(0)))

        for nref=0:maxref
            res=run_fet(nref=nref,conv=true)

            for case in cases
                idrain=res[case]["Idrain"]
                diff=norm(refcurr-idrain,2)
                err=alldata[case]["error"]
                push!(err,diff)
            end
        end
        save("conv.jld2",alldata)
        return
    else
        alldata=Dict()
        alldata=load("conv.jld2")
    end

    PyPlot.clf()
    X=nothing
    Y=nothing
    Y1=nothing
    for case in cases
        err=alldata[case]["error"]
        L=0.1
        X=[L/(10*2^i) for i=1:length(err)]
        Y=[150*X[i] for i=1:length(err)]
        line=alldata[case]["line"]
        flux=alldata[case]["name"]
        marker=alldata[case]["marker"]
        markersize=alldata[case]["markersize"]
        PyPlot.loglog(X,err,line,marker="o",label=flux)
    end
    PyPlot.xlim(5.0e-5,1.0e-2)
    PyPlot.loglog(X,Y,"k--",label="\$O(h)\$")
    PyPlot.loglog(X,100*Y.^2,"k-.",label="\$O(h^2)\$")
    PyPlot.ylabel("\$||\\mathrm{IV}-\\mathrm{IV}_{ref}||_2\$")
    PyPlot.xlabel("\$h_x\$")
    PyPlot.grid()
    PyPlot.legend(loc="lower right")

    show()

    mysavefig("fet-conv.pdf")
end
```

# Helper functions

## Create all figures in the paper

By default, parameter `quick` is true, resulting
in ca 5min of calculations on a laptop. This leads to
to convergence results in a smaller range of
refinement levels. With quick=true, calculations run
overnight (mostly due to Fig.7).

```julia
function makefigures(;quick=true)
    for i=1:7
        makefigure(i,quick=quick)
    end
end
```

## Calculations for given figure

```julia
function makefigure(ifig; quick=true)

    n1d=100
    if quick
        n1d=75
    end
    ifig==1 && run_to_equilibrium(n=n1d,subcase="1",genplot=true, cases=["c","a","b","s"])
    ifig==2 && run_to_equilibrium(n=n1d,subcase="2",genplot=true, cases=["c","a","b","s"])
    ifig==3 && run_to_equilibrium(n=n1d,subcase="3",genplot=true, cases=["c","a","b","s"])


    maxref=12
    if quick
        maxref=8
    end
    ifig==4 && run_compare(maxref=maxref,stationary=false)
    ifig==5 && run_compare(maxref=maxref,stationary=true)

    nref=3
    if quick
        nref=1
    end
    ifig==6 && run_fet(nref=nref, genplot=true)

    if ifig==7
        nref=7
        maxref=6
        if (quick)
            nref=4
            maxref=3
        end
```

       run_fet(nref=nref, genref=true)
       run_conv(maxref=maxref,genconv=true)

```julia
        run_conv(maxref=maxref,genconv=false)
    end
end
```

## Wrap savefig in order to print out name

```julia
function mysavefig(figname)
    PyPlot.savefig(figname)
    println("Saved figure $(figname)")
end
```

## Prepare dictionary for collecting results

```julia
function DataDict()
    alldata=Dict()

    push!(alldata,Pair("c",Dict()))
    push!(alldata["c"],Pair("name","Centered"))
    push!(alldata["c"],Pair("flux",cflux!))
    push!(alldata["c"],Pair("line","m-"))
    push!(alldata["c"],Pair("marker","o"))
    push!(alldata["c"],Pair("markersize",4))


    push!(alldata,Pair("a",Dict()))
    push!(alldata["a"],Pair("name","Activity"))
    push!(alldata["a"],Pair("flux",aflux!))
    push!(alldata["a"],Pair("line","g-"))
    push!(alldata["a"],Pair("marker","x"))
    push!(alldata["a"],Pair("markersize",10))


    push!(alldata,Pair("b",Dict()))
    push!(alldata["b"],Pair("name","Bess-Ch"))
    push!(alldata["b"],Pair("flux",bflux!))
    push!(alldata["b"],Pair("line","r-"))
    push!(alldata["b"],Pair("marker","+"))
    push!(alldata["b"],Pair("markersize",10))





    push!(alldata,Pair("s",Dict()))
    push!(alldata["s"],Pair("name","Sedan"))
    push!(alldata["s"],Pair("flux",sflux!))
    push!(alldata["s"],Pair("line","b-"))
    push!(alldata["s"],Pair("marker","o"))
    push!(alldata["s"],Pair("markersize",4))

    return alldata
end
```

## Plot 1D solution

```julia
function plot_solution(sys,U0)
    ildata=data(sys)
    iΦ=ildata.iΦ
    ic=ildata.ic
    PyPlot.clf()
    @views begin
        PyPlot.plot(sys.grid.coord[1,:],U0[iΦ,:], label="\$\\Phi\$", color="g")
        PyPlot.plot(sys.grid.coord[1,:],U0[ic,:], label="c", color="b")
    end
```

   PyPlot.ylim(0,2)

```julia
    PyPlot.xlabel("\$x\$")
    PyPlot.grid()
    PyPlot.legend(loc="upper left")
    PyPlot.pause(1.0e-10)
end
```

## Plot 2D solution

```julia
function plot2d(data, X,Y,sol;umax=25,ysc="linear", xlbl="x",ylbl="y", uraster=1, casename="")
    iΦ=data.iΦ
    ic=data.ic
    PyPlot.clf()
    umax=uraster*ceil(umax/uraster)
    v=collect(-umax:umax/100:umax)
    t=collect(-umax:umax/5:umax)

    subplot(121)
    @views title(@sprintf("\$\\Phi\$: min=%.4f max=%.4f",minimum(sol[iΦ,:]),maximum(sol[iΦ,:])))
    cnt=contourf(X,Y,transpose(reshape(sol[iΦ,:],length(X),length(Y))),v,cmap=ColorMap("bwr"), antialised=false)
    for c in cnt.collections
        c.set_edgecolor("face")
    end

    colorbar(ticks=t,boundaries=v)
    contour(X,Y,transpose(reshape(sol[iΦ,:],length(X),length(Y))),colors="k",linetypes="-",t,linewidths=0.75)

    yscale(ysc)
    xlabel(xlbl)
    ylabel(ylbl)

    subplot(122)
    v=collect(-0.05:0.01:1.05)
    t=collect(0:0.1:1)

    mmax=maximum(sol[ic,:])
    if mmax>0.99
        @views title(@sprintf("c: min=%.3e max=1-%.3e",minimum(sol[ic,:]),1-mmax))
    else
        @views title(@sprintf("c: min=%.3e max=%.3e",minimum(sol[ic,:]),mmax))
    end


    cnt=contourf(X,Y,transpose(reshape(sol[ic,:],length(X),length(Y))),v, cmap=ColorMap("hot") )
    colorbar(ticks=t,boundaries=v)
    for c in cnt.collections
        c.set_edgecolor("face")
    end
    contour(X,Y,transpose(reshape(sol[ic,:],length(X),length(Y))),colors="k",linetypes="-",t,linewidths=0.75)
    yscale(ysc)
    xlabel(xlbl)
    ylabel(ylbl)
    tight_layout()
    pause(1.0e-10)
    show()
end
```

## Make pdf from this source code

```julia
function makepdf()
    Literate.markdown("UnipolarDriftDiffusion.jl", ".",documenter=false)
    run(`pandoc UnipolarDriftDiffusion.md -t latex --metadata-file=header.yaml --number-sections --toc  --highlight-style tango --pdf-engine=xelatex -o UnipolarDriftDiffusion.pdf`)
end
```

## Generate readme.md from this source code

```julia
function makereadme()
    Literate.markdown("UnipolarDriftDiffusion.jl", ".",documenter=false)
    run(`mv UnipolarDriftDiffusion.md README.md`)
end
```

# Calculations not used in the paper

## Double layer capacitance

```julia
function run_dlcap(;physics=DefaultPhysics(),n=100,pyplot=false)
    iΦ=physics.data.iΦ
    ic=physics.data.ic
    L=25
    h=L/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:L))

    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])

    sys.boundary_values[iΦ,1]=0
    sys.boundary_values[iΦ,2]=0.0

    sys.boundary_factors[iΦ,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[iΦ,2]=VoronoiFVM.Dirichlet

    sys.boundary_values[ic,2]=0.5
    sys.boundary_factors[ic,2]=VoronoiFVM.Dirichlet


    control=VoronoiFVM.NewtonControl()
    control.verbose=true
    control.damp_initial=0.001
    print("calculating double layer capacitance")
    delta=1.0e-4

    inival=unknowns(sys)
    @views inival[iΦ,:].=0
    @views inival[ic,:].=0.5

    sys.boundary_values[iΦ,1]=0

    dphi=1.0e-1
    phimax=10
    delta=1.0e-4
    vplus=zeros(0)
    cdlplus=zeros(0)
    vminus=zeros(0)
    cdlminus=zeros(0)
    if pyplot
        plot_solution(sys,inival)
    end
    inival0=copy(inival)
    sol=copy(inival)
    for dir in [1,-1]
        inival.=inival0
        phi=0.0
        while phi<phimax
            sys.boundary_values[iΦ,1]=dir*phi
            solve!(sol,inival,sys)
            Q=integrate(sys,physics.reaction,sol)
            sys.boundary_values[iΦ,1]=dir*phi+delta
            inival.=sol

            solve!(sol,inival,sys)
            if pyplot
                plot_solution(sys,sol)
            end
            Qdelta=integrate(sys,physics.reaction,sol)
            cdl=(Qdelta[iΦ]-Q[iΦ])/delta
            if dir==1
                push!(vplus,dir*phi)
                push!(cdlplus,cdl)
            else
                push!(vminus,dir*phi)
                push!(cdlminus,cdl)
            end
            phi+=dphi
        end
    end
    if pyplot
        PyPlot.clf()
        PyPlot.plot(vplus,cdlplus,color="g")
        PyPlot.plot(vminus,cdlminus,color="g")
        PyPlot.grid()
        PyPlot.legend(loc="upper right")
        PyPlot.pause(1.0e-10)
    end
end
```

## Simple time evolution

```julia
function run_transient(;L=50,physics=DefaultPhysics(),n=100,pyplot=false)
    h=L/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:L))
    data=physics.data
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])

    sys.boundary_values[data.iΦ,1]=0.0
    sys.boundary_values[data.iΦ,2]=0.0


    sys.boundary_factors[data.iΦ,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[data.iΦ,2]=VoronoiFVM.Dirichlet

    sys.boundary_values[data.ic,2]=0.99
    sys.boundary_values[data.ic,1]=0.01
    sys.boundary_factors[data.ic,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[data.ic,2]=VoronoiFVM.Dirichlet

    inival=unknowns(sys)
    @views inival[data.iΦ,:].=0
    @views inival[data.ic,:].=0.5
    U=unknowns(sys)

    control=VoronoiFVM.NewtonControl()
    control.verbose=true
    print("time loop")
    control.max_lureuse=0
    control.damp_initial=0.01
    t=0.0
    tend=1.0
    tstep=1.0e-4
    while t<tend
        t=t+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        @printf("time=%g\n",t)
        if pyplot
            plot_solution(sys,U)
        end
        tstep*=1.4
    end
end

end
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

