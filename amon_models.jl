using CSV
using DataFrames
using Missings
using Gadfly

#Code to retrieve dataset from Pangaea
function get_pangaea_data(url)
    lines = readlines(download(url))
    goodlines = []
    commentline = false
    for line in lines
        if startswith(line, "/*")
            commentline = true
            continue
        elseif commentline
            if startswith(line, "*/")
                commentline = false
                continue
            else
                continue
            end
        elseif commentline == false
            push!(goodlines, line)
        end
    end
    colnames = split(goodlines[1], "\t")
    data = permutedims(stack(map(l -> split(l, "\t"),goodlines[2:end])),(2,1))
    df = DataFrame(data, :auto)
    rename!(df, colnames)
    @show df
    return df
end

#Download data
med = get_pangaea_data("https://doi.pangaea.de/10.1594/PANGAEA.809994?format=textfile")
med = ifelse.(med .== "", missing, med)
med[!,1:5] = passmissing(parse).(Float64, med[!,1:5])

#Get measurement O2 and pH
med = med[ismissing.(med.pH) .== 0 .&& ismissing.(med[!,"O2 [µmol/l]"]) .== 0,:]

#Get values from depth 0 and down
med = med[med[!,Symbol("Depth sed [m]")] .>0.0,:]

#Convert all units to molar
med[!,:O2_mol] = med[!,Symbol("O2 [µmol/l]")] ./ 1E6
med[!,:H_mol] = exp10.(-med[!,:pH])

#Convert 0 oxygen to 1 nanomolar (1E-9 M)
med[med[!,:O2_mol] .== 0.0,:O2_mol] .= 1E-9

function deltaG(T = Float64{}, Q = Float64{})
    #output delta G in kJ
    R = 8.314462618
    return R*T*log(ℯ, Q) / 1000
end

F = 96485.3321233100184
voltage_drop_mm = -0.013 #V mm-1
z = 4 #4 electrons transfered per O2

Gdrop_mm = -z*F*voltage_drop_mm/1000 #delta G drop per mm in kJ

#get sediment depth in mm
med[!, :Depth_mm] = med[!,Symbol("Depth sed [m]")].*1000

#calculate reaction quotient
med[!, :Q] =  (med[!, :O2_mol]*med[1, :H_mol]^4)./(med[1, :O2_mol]*med[!, :H_mol].^4)

#calculate delta G
med[!, :deltaG] = deltaG.(293, med[!, :Q])
med[!, :loss] = med[!, :Depth_mm] .* Gdrop_mm
med[!, :netG] = med[!, :deltaG] .+ med[!, :loss]


@show med

O2_conc_plot = plot(med, x = Symbol("O2 [µmol/l]"), y = :Depth_mm, 
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Guide.ylabel("Depth (mm)"),
    Guide.xlabel("O2 (µM)")
    )

pH_plot = plot(med, x = :pH, y = :Depth_mm, 
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Guide.ylabel(""),
    Guide.xlabel("pH")
    )

deltaG_plot = plot(med, x = :netG, y = :Depth_mm, xintercept = [-10,0.0],
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Geom.vline(color = ["red","black"]),
    Guide.xlabel("ΔG (kJ/mol)"),
    Guide.ylabel("")
    )

hstack(O2_conc_plot, pH_plot, deltaG_plot)