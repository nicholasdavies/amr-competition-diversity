function all(distribution)
    return string.rep(distribution, 27, ",")
end

function sweep(parameter, pmin, pmax, trials)
    return "trial -> '" .. parameter .. "', -1, " .. pmin .. " + trial * " .. (pmax - pmin) .. " / " .. (trials - 1)
end

function country_allocate(param_name, ...)
    x = table.pack(...)
    ret = {}
    for i = 0, 26 do
        for j = 0, 9 do
            table.insert(ret, param_name);
            table.insert(ret, i * 10 + j);
            table.insert(ret, x[i + 1]);
        end
     end
     return table.unpack(ret)
end

prior_beta = "G 5 0.35 I 1.42 1.42"
prior_u = "G 5 0.1625 I 0.65 0.65"
prior_c = "B 1.5 8.5 I 0.20 0.20"
prior_g = "B 10.85 1.15 I 0.985 0.985 T 0 0.999999"
prior_g_within = "B 10.85 1.15 I 0.985 0.985 T 0 0.999999"
prior_b = "G 2 0.125 I 0.19 0.19"
prior_k = "N 1 0.5 T 0 50 I 1.7 1.7"
prior_a = "G 2 5 I 20 20 T 0 100"
prior_delta = "B 20 25 I 0.5 0.5"
prior_shape = "G 4 2 I 2.5 2.5"
prior_z = "L 0.16 0.4 I 1 1"