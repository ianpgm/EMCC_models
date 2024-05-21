using DataFrames
using Gadfly, ColorSchemes

function twopart_concentration_profile(
    upperconc = Float64{},
    lowerconc = Float64{},
    penetration_distance = Float64{},
    distance = Float64{})
    if distance <= penetration_distance
        conc = upperconc - (distance / penetration_distance)*(upperconc - lowerconc)
    elseif distance > penetration_distance
        conc = lowerconc
    end
    return conc
end

function deltaG(T = Float64{}, Q = Float64{})
    #output delta G in kJ
    R = 8.314462618
    return R*T*log(ℯ, Q) / 1000
end

function Gdrop_mm(z)
    #z is number of electrons transferred
    F = 96485.3321233100184
    voltage_drop_mm = -0.013 #V mm-1

    Gdrop_mm = -z*F*voltage_drop_mm/1000 #delta G drop per mm in kJ

    return Gdrop_mm
end

function make_model(z, penetration_distance, conc_diff_orders)
    lowerconc = 1E-9
    upperconc_exponent = log(10,lowerconc)+conc_diff_orders
    upperconc = 10^upperconc_exponent

    model = DataFrame(Depth = 0:0.1:20)
    model[!, :acceptor] = twopart_concentration_profile.(upperconc, lowerconc, penetration_distance, model[!, :Depth])
    model[!, :Q] =  model[!, :acceptor] ./ upperconc
    model[!, :deltaG] = deltaG.(293, model[!, :Q])
    model[!, :loss] = (model[!, :Depth]) .* Gdrop_mm(z)
    model[!, :netG] = model[!, :deltaG] .+ model[!, :loss]
    return model
end

function min_delta_G(model)
    return minimum(model[!,:netG])
end

function mm_below_threshold(model, threshold)
    return count(model[!,:netG] .< threshold)/10
end

O2_model = make_model(4, 1, 9)

mm_below_threshold(O2_model, -10)
min_delta_G(O2_model)

model_df = DataFrame(electrons = [], orders_magnitude = [], minimum_deltaG = [], mm_below_threshold = [])

for e in 2:2:8
    for o in 0:9
        newmodel = make_model(e, 1, o)
        minimum_deltaG = min_delta_G(newmodel)
        mm_under_threshold = mm_below_threshold(newmodel, -10)
        newline = DataFrame(electrons = e, orders_magnitude = o, minimum_deltaG = minimum_deltaG, mm_below_threshold = mm_under_threshold)
        model_df = vcat(model_df, newline)
    end
end

G = plot(model_df, x = :orders_magnitude, y = :electrons, colour = :minimum_deltaG,
    Geom.rectbin(),
    Coord.cartesian(xmin=-0.5, xmax=9.5),
    Guide.ylabel("Electrons transferred"),
    Guide.xlabel("Log concentration difference"),
    Guide.colorkey(""),
    Theme(grid_line_width=0mm),
    Guide.yticks(label=true, ticks=[2:2:8...], orientation=:horizontal),
    Guide.xticks(label=true, ticks=[0:9...], orientation=:horizontal),
    Scale.color_continuous(colormap = viridis_cm_rev)
    )

L = plot(model_df, x = :orders_magnitude, y = :electrons, colour = :mm_below_threshold,
    Geom.rectbin(),
    Coord.cartesian(xmin=-0.5, xmax=9.5),
    Guide.ylabel("Electrons transferred"),
    Guide.xlabel("Log concentration difference"),
    Guide.colorkey(""),
    Theme(grid_line_width=0mm),
    Guide.yticks(label=true, ticks=[2:2:8...], orientation=:horizontal),
    Guide.xticks(label=true, ticks=[0:9...], orientation=:horizontal),
    Scale.color_continuous(minvalue = 0, maxvalue = 16, colormap = viridis_cm)
    )

vstack(G,L)


example_model = make_model(4, 1, 6)

@show example_model

conc_plot = plot(example_model, x = :acceptor, y = :Depth, 
    Geom.path(),
    Coord.Cartesian(yflip = true, ymax = 8),
    Guide.ylabel("Depth (mm)"),
    Guide.xlabel("Electron acceptor (M)")
    )

deltaG_plot = plot(example_model, x = :netG, y = :Depth, xintercept = [-10.0, 0],
    Geom.path(),
    Coord.Cartesian(yflip = true, ymax = 8, xmax = 10),
    Geom.vline(color = "black"),
    Guide.xlabel("ΔG (kJ/mol)"),
    Guide.ylabel("")
    )

hstack(conc_plot, deltaG_plot)


