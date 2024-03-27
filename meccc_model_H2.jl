using DataFrames
using Gadfly

function twopart_concentration_profile(
    upperconc = Float64{},
    lowerconc = Float64{},
    penetration_distance = Float64{},
    top_distance = Float64{},
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
z = 2 #2 electrons transfered per H2

Gdrop_mm = -z*F*voltage_drop_mm/1000 #delta G drop per mm in kJ

H2_model = DataFrame(Depth = 0:0.01:5)
H2_model[!, :H2] = twopart_concentration_profile.(1E-3, 1E-9, 1, 5, H2_model[!, :Depth])
H2_model[!, :H2_mM] = H2_model[!, :H2] * 1000
H2_model[!, :Q] = H2_model[!, :H2] ./ 1E-3
H2_model[!, :deltaG] = deltaG.(293, H2_model[!, :Q])
H2_model[!, :loss] = (H2_model[!, :Depth]) .* Gdrop_mm
H2_model[!, :netG] = H2_model[!, :deltaG] .+ H2_model[!, :loss]

conc_plot = plot(H2_model, x = :H2_mM, y = :Depth, 
    Geom.path(),
    #Coord.Cartesian(yflip = true),
    Guide.ylabel("Distance from H2 source (mm)"),
    Guide.xlabel("H2 (mM)")
    )

deltaG_plot = plot(H2_model, x = :netG, y = :Depth, xintercept = [0.0],
    Geom.path(),
    #Coord.Cartesian(yflip = true),
    Geom.vline(color = "black"),
    Guide.xlabel("ΔG (kJ/mol)"),
    Guide.ylabel("")
    )

hstack(conc_plot, deltaG_plot)


#= H2_model_stacked = stack(H2_model, [:H2, :netG])

p = plot(H2_model_stacked, 
    x = :value, y = :Depth, xgroup = :variable,
    Geom.subplot_grid(Geom.path, free_x_axis = true),
    Guide.ylabel("Distance from H2 source (mm)"),
    Guide.xlabel("")
    )

p |> SVG("H2.svg", 6inch, 4inch) =#