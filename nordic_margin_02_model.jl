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
nordicmargin = get_pangaea_data("https://doi.pangaea.de/10.1594/PANGAEA.772691?format=textfile")
nordicmargin = ifelse.(nordicmargin .== "", missing, nordicmargin)
nordicmargin[!,1:5] = passmissing(parse).(Float64, nordicmargin[!,1:5])

#Get measurement O2 and pH
nordicmargin = nordicmargin[ismissing.(nordicmargin.pH) .== 0 .&& ismissing.(nordicmargin[!,"O2 [µmol/l]"]) .== 0,:]

#Get values from depth 0 and down
nordicmargin = nordicmargin[nordicmargin[!,Symbol("Depth sed [m]")] .>0.0,:]

#Convert all units to molar
nordicmargin[!,:O2_mol] = nordicmargin[!,Symbol("O2 [µmol/l]")] ./ 1E6
nordicmargin[!,:H_mol] = exp10.(-nordicmargin[!,:pH])

#Convert oxygen at deeper than 2.7 cm to 1 nanomolar (1E-9 M)
nordicmargin[nordicmargin[!,Symbol("Depth sed [m]")] .> 0.02775,:O2_mol] .= 1E-9

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
nordicmargin[!, :Depth_mm] = nordicmargin[!,Symbol("Depth sed [m]")].*1000

#calculate reaction quotient
nordicmargin[!, :Q] =  (nordicmargin[!, :O2_mol]*nordicmargin[1, :H_mol]^4)./(nordicmargin[1, :O2_mol]*nordicmargin[!, :H_mol].^4)

#calculate delta G
nordicmargin[!, :deltaG] = deltaG.(293, nordicmargin[!, :Q])
nordicmargin[!, :loss] = nordicmargin[!, :Depth_mm] .* Gdrop_mm
nordicmargin[!, :netG] = nordicmargin[!, :deltaG] .+ nordicmargin[!, :loss]


@show nordicmargin

O2_conc_plot = plot(nordicmargin, x = Symbol("O2 [µmol/l]"), y = :Depth_mm, 
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Guide.ylabel("Depth (mm)"),
    Guide.xlabel("O2 (µM)")
    )

pH_plot = plot(nordicmargin, x = :pH, y = :Depth_mm, 
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Guide.ylabel(""),
    Guide.xlabel("pH")
    )

deltaG_plot = plot(nordicmargin, x = :netG, y = :Depth_mm, xintercept = [-10,0.0],
    Geom.path(),
    Coord.Cartesian(yflip = true),
    Geom.vline(color = ["red","black"]),
    Guide.xlabel("ΔG (kJ/mol)"),
    Guide.ylabel("")
    )

hstack(O2_conc_plot, pH_plot, deltaG_plot)