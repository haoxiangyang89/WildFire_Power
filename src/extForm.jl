# model and solve the extensive formulation

function extForm(netData,scenData,sol_ind = false)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)));
    
    if sol_ind
        optimize!(mp);
        return objective_value(mp);
    else
        return mp;
    end
end
