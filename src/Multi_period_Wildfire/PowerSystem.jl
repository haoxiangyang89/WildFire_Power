using JuMP
import DataFrames
import GLPK
import Plots
import StatsPlots



## prepare data

function ThermalGenerator(
    min::Float64,
    max::Float64,
    fixed_cost::Float64,
    variable_cost::Float64,
)
    return (
        min = min,
        max = max,
        fixed_cost = fixed_cost,
        variable_cost = variable_cost,
    )
end

generators = [
    ThermalGenerator(0.0, 1000.0, 1000.0, 50.0),
    ThermalGenerator(300.0, 1000.0, 0.0, 100.0),
]


WindGenerator(variable_cost::Float64) = (variable_cost = variable_cost,)

wind_generator = WindGenerator(50.0)


function Scenario(demand::Float64, wind::Float64)
    return (demand = demand, wind = wind)
end

scenario = Scenario(1500.0, 200.0)




## build model
function solve_ed(generators::Vector, wind, scenario)
    # Define the economic dispatch (ED) model
    ed = Model(GLPK.Optimizer)
    # Define decision variables
    # power output of generators
    N = length(generators)
    @variable(ed, generators[i].min <= g[i = 1:N] <= generators[i].max)
    # wind power injection
    @variable(ed, 0 <= w <= scenario.wind)
    # Define the objective function
    @objective(
        ed,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w,
    )
    # Define the power balance constraint
    @constraint(ed, sum(g[i] for i in 1:N) + w == scenario.demand)
    # Solve statement
    optimize!(ed)
    # return the optimal value of the objective function and its minimizers
    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(ed),
    )
end




## solve the model
solution = solve_ed(generators, wind_generator, scenario);

println("Dispatch of Generators: ", solution.g, " MW")
println("Dispatch of Wind: ", solution.w, " MW")
println("Wind spillage: ", solution.wind_spill, " MW")
println("Total cost: \$", solution.total_cost)




#############################################################################################################################
############################### Economic dispatch with adjustable incremental costs #########################################
#############################################################################################################################

function scale_generator_cost(g, scale)
    return ThermalGenerator(g.min, g.max, g.fixed_cost, scale * g.variable_cost)
end

start = time()
c_g_scale_df = DataFrames.DataFrame(
    # Scale factor
    scale = Float64[],
    # Dispatch of Generator 1 [MW]
    dispatch_G1 = Float64[],
    # Dispatch of Generator 2 [MW]
    dispatch_G2 = Float64[],
    # Dispatch of Wind [MW]
    dispatch_wind = Float64[],
    # Spillage of Wind [MW]
    spillage_wind = Float64[],
    # Total cost [$]
    total_cost = Float64[],
)

for c_g1_scale in 0.5:0.1:3.0
    # Update the incremental cost of the first generator at every iteration.
    new_generators = scale_generator_cost.(generators, [c_g1_scale, 1.0])
    # Solve the ed problem with the updated incremental cost
    sol = solve_ed(new_generators, wind_generator, scenario)
    push!(
        c_g_scale_df,
        (c_g1_scale, sol.g[1], sol.g[2], sol.w, sol.wind_spill, sol.total_cost),
    )
end
print(string("elapsed time: ", time() - start, " seconds"))





#############################################################################################################################
########################################## Modifying the JuMP model in-place ################################################
#############################################################################################################################


function solve_ed_inplace(
    generators::Vector,
    wind,
    scenario,
    scale::AbstractVector{Float64},
)
    obj_out = Float64[]
    w_out = Float64[]
    g1_out = Float64[]
    g2_out = Float64[]
    # This function only works for two generators
    @assert length(generators) == 2
    ed = Model(GLPK.Optimizer)
    set_silent(ed)
    N = length(generators)
    @variable(ed, generators[i].min <= g[i = 1:N] <= generators[i].max)
    @variable(ed, 0 <= w <= scenario.wind)
    @objective(
        ed,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w,
    )
    @constraint(ed, sum(g[i] for i in 1:N) + w == scenario.demand)
    for c_g1_scale in scale
        @objective(
            ed,
            Min,
            c_g1_scale * generators[1].variable_cost * g[1] +
            generators[2].variable_cost * g[2] +
            wind.variable_cost * w,
        )
        optimize!(ed)
        push!(obj_out, objective_value(ed))
        push!(w_out, value(w))
        push!(g1_out, value(g[1]))
        push!(g2_out, value(g[2]))
    end
    df = DataFrames.DataFrame(
        scale = scale,
        dispatch_G1 = g1_out,
        dispatch_G2 = g2_out,
        dispatch_wind = w_out,
        spillage_wind = scenario.wind .- w_out,
        total_cost = obj_out,
    )
    return df
end

start = time()
inplace_df = solve_ed_inplace(generators, wind_generator, scenario, 0.5:0.1:3.0)
print(string("elapsed time: ", time() - start, " seconds"))





#############################################################################################################################
######################################### Inefficient usage of wind generators ##############################################
#############################################################################################################################

demand_scale_df = DataFrames.DataFrame(
    demand = Float64[],
    dispatch_G1 = Float64[],
    dispatch_G2 = Float64[],
    dispatch_wind = Float64[],
    spillage_wind = Float64[],
    total_cost = Float64[],
)

function scale_demand(scenario, scale)
    return Scenario(scale * scenario.demand, scenario.wind)
end

for demand_scale in 0.2:0.1:1.4
    new_scenario = scale_demand(scenario, demand_scale)
    sol = solve_ed(generators, wind_generator, new_scenario)
    push!(
        demand_scale_df,
        (
            new_scenario.demand,
            sol.g[1],
            sol.g[2],
            sol.w,
            sol.wind_spill,
            sol.total_cost,
        ),
    )
end

demand_scale_df




dispatch_plot = StatsPlots.@df(
    demand_scale_df,
    Plots.plot(
        :demand,
        [:dispatch_G1, :dispatch_G2],
        labels = ["G1" "G2"],
        title = "Thermal Dispatch",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand",
        ylabel = "Dispatch [MW]",
    ),
)

wind_plot = StatsPlots.@df(
    demand_scale_df,
    Plots.plot(
        :demand,
        [:dispatch_wind, :spillage_wind],
        labels = ["Dispatch" "Spillage"],
        title = "Wind",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand [MW]",
        ylabel = "Energy [MW]",
    ),
)

Plots.plot(dispatch_plot, wind_plot)












#############################################################################################################################
###################################################### Unit Commitment ######################################################
#############################################################################################################################

function solve_uc(generators::Vector, wind, scenario)
    uc = Model(GLPK.Optimizer)
    set_silent(uc)
    N = length(generators)
    @variable(uc, generators[i].min <= g[i = 1:N] <= generators[i].max)
    @variable(uc, 0 <= w <= scenario.wind)
    @constraint(uc, sum(g[i] for i in 1:N) + w == scenario.demand)
    # !!! New: add binary on-off variables for each generator
    @variable(uc, u[i = 1:N], Bin)
    @constraint(uc, [i = 1:N], g[i] <= generators[i].max * u[i])
    @constraint(uc, [i = 1:N], g[i] >= generators[i].min * u[i])
    @objective(
        uc,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w +
        # !!! new
        sum(generators[i].fixed_cost * u[i] for i in 1:N)
    )
    optimize!(uc)
    status = termination_status(uc)
    if status != OPTIMAL
        return (status = status,)
    end
    return (
        status = status,
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        u = value.(u),
        total_cost = objective_value(uc),
    )
end



solution = solve_uc(generators, wind_generator, scenario)

println("Dispatch of Generators: ", solution.g, " MW")
println("Commitments of Generators: ", solution.u)
println("Dispatch of Wind: ", solution.w, " MW")
println("Wind spillage: ", solution.wind_spill, " MW")
println("Total cost: \$", solution.total_cost)