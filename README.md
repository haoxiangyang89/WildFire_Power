# Power System Risk Minimization under Wildfire Disruptions

A Julia code that simulates the wildfire scenarios and optimizes the power system risk minimization problem. 

## 

### Data 

The required data for the problem is as follows: 

- $\texttt{indexSets}$: the network index set, including the set of buses $\mathcal{B}$ and set of lines $\mathcal{L}$, etc.;
- $\texttt{paramOPF}$: the necessary parameters for optimal power flow problem;
- $\texttt{paramDemand}$: the necessary parameters for the load demand;
- $\Omega\_\texttt{rv}$ and $\texttt{prob}$: the set of scenarios and the probability of each scenario $\omega \in \Omega\_\texttt{rv}$.

You can generate data from the file in $\texttt{data}$:

```julia
.../WildFire_Power/data/testData_RTS/generateData_RTS_GMLC.jl 
```



### Decomposition algorithm:

Go to the file:

```julia
.../src/alg/loadMod.jl
```

Upload all the required Julia packages and files. Then you can call the function: 

```julia
sddipResult = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
```

