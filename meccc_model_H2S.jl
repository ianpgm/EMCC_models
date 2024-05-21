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
z = 8 #8 electrons transfered per H2S

Gdrop_mm = -z*F*voltage_drop_mm/1000 #delta G drop per mm in kJ

#H2S in sediment 5 µM

H2S_model = DataFrame(Depth = 0:0.01:10)
H2S_model[!, :H2S] = twopart_concentration_profile.(1E-9, 2E-3, 7, H2S_model[!, :Depth])
H2S_model[!, :H2S_mM] = H2S_model[!, :H2S] .* 1E3
H2S_model[!, :Q] =  1E-9 ./ H2S_model[!, :H2S]
H2S_model[!, :deltaG] = deltaG.(293, H2S_model[!, :Q])
H2S_model[!, :loss] = H2S_model[!, :Depth] .* Gdrop_mm
H2S_model[!, :netG] = H2S_model[!, :deltaG] .+ H2S_model[!, :loss]

#H2S_model_stacked = stack(H2S_model, [:H2S, :netG])

conc_plot = plot(H2S_model, x = :H2S_mM, y = :Depth, 
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Guide.ylabel("Depth (mm)"),
    Guide.xlabel("H2S (mM)")
    )

deltaG_plot = plot(H2S_model, x = :netG, y = :Depth, xintercept = [0.0],
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Geom.vline(color = "black"),
    Guide.xlabel("ΔG (kJ/mol)"),
    Guide.ylabel("")
    )

hstack(conc_plot, deltaG_plot)


#plot(H2S_model_stacked, 
#    x = :value, y = :Depth, xgroup = :variable,
#    Geom.subplot_grid(Geom.path, free_x_axis = true, yflip=true),
#    Guide.ylabel("Depth (mm)"),
#    Guide.xlabel("")
#    )

#p |> SVG("O2_01mm.svg", 6inch, 4inch)