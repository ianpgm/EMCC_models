using DataFrames
using Gadfly

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

F = 96485.3321233100184
voltage_drop_mm = -0.013 #V mm-1
z = 4 #4 electrons transfered per O2

Gdrop_mm = -z*F*voltage_drop_mm/1000 #delta G drop per mm in kJ

#O2 max. solubility 0.00026 M 

O2_model = DataFrame(Depth = 0:0.1:5)
O2_model[!, :O2] = twopart_concentration_profile.(0.0003, 1E-6, 0.5, O2_model[!, :Depth])
O2_model[!, :O2_micromol] =  O2_model[!, :O2] *1E6
O2_model[!, :Q] =  O2_model[!, :O2] ./ 0.0003
O2_model[!, :deltaG] = deltaG.(293, O2_model[!, :Q])
O2_model[!, :loss] = (O2_model[!, :Depth]) .* Gdrop_mm
O2_model[!, :netG] = O2_model[!, :deltaG] .+ O2_model[!, :loss]


conc_plot = plot(O2_model, x = :O2_micromol, y = :Depth, 
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Guide.ylabel("Depth (mm)"),
    Guide.xlabel("O2 (µM)")
    )

deltaG_plot = plot(O2_model, x = :netG, y = :Depth, xintercept = [0.0],
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Geom.vline(color = "black"),
    Guide.xlabel("ΔG (kJ/mol)"),
    Guide.ylabel("")
    )

hstack(conc_plot, deltaG_plot)

#= 

O2_model_stacked = stack(O2_model, [:O2, :netG])

p = plot(O2_model_stacked, 
    x = :value, y = :Depth, xgroup = :variable,
    Geom.subplot_grid(Geom.path, free_x_axis = true),
    Guide.ylabel("Distance from base of model (mm)"),
    Guide.xlabel("")
    )

p |> SVG("O2_05mm.svg", 6inch, 4inch)

=#