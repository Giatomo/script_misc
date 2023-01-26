using Gadfly: ismissing
##

using DataFrames
using CSV
using ShiftedArrays
using Statistics
using Gadfly
using IntervalArithmetic
using Chain
using Base: argmin
##
@inline function allequal(x)
    length(x) < 2 && return true, x[1]
    e1 = x[1]
    @inbounds for i=2:length(x)
        x[i] == e1 || return false, nothing
    end
    return true, x[1]
end


function radialcoordinate(Δx::Number, Δy::Number)
    return √(Δx^2 + Δy^2)
end

function radialcoordinate(Δx::Missing, Δy::Missing)
    return(missing)
end

function angularcoordinate(Δx::Number, Δy::Number)
    θ = atan(Δy/Δx) * (180/π)
    quadrant = getquadrant(Δy, Δx)
    # Quadrant I => Do nothing
    if quadrant == 1
        return θ

    # Quadrant IV => Add 360°
    elseif quadrant == 4  
        return θ + 360

    # Quadrant II/III => Add 180°
    else  
        return θ + 180
    end
end


function getquadrant(Δx::Number, Δy::Number)
    if (Δx >= 0 && Δy >= 0) # Quadrant I
        return 1
    elseif (Δx < 0 && Δy >= 0) # Quadrant II
        return 2
    elseif (Δx < 0 && Δy < 0) # Quadrant III
        return 3
    else  # Quadrant IV 
        return 4
    end
end

function angularcoordinate(Δx::Missing, Δy::Missing)
    return missing
end

function isoppositedirection(α::Number, θ::Number, interval::Number)
    return abs(α - θ) ∈ 180±interval ? -1 : 1
end

##
bacteria = DataFrame(CSV.File("/home/thomas/Bureau/BACT.csv"))
spots = DataFrame(CSV.File("/home/thomas/Bureau/MAX_1.csv"))

##
tidy_bacteria = @chain bacteria begin
    rename(
        :NAME => :B_ID,
        :LOCATION_X => :B_X,
        :LOCATION_Y => :B_Y, 
        :TRAJECTORY => :T_ID,
        :POSITION => :FRAME,
        :TRAJECTORY_START => :FRAME_START,
        :TRAJECTORY_END => :FRAME_END, 
        :MAXIMA => :N_MAXIMA)
    select(
        :B_ID,
        :B_X, 
        :B_Y,
        :T_ID,
        :FRAME, 
        :FRAME_START, 
        :FRAME_END, 
        :N_MAXIMA)
end
##

tidy_spot = @chain spots begin
    rename(
        :NAME => :M_ID,
        :LOCATION_X => :M_X,
        :LOCATION_Y => :M_Y,
        :PARENT => :B_ID)
    select(
        :B_ID, 
        :M_ID, 
        :M_X, 
        :M_Y, 
        :INTENSITY, 
        :ZSCORE)
end
##
df = @chain tidy_bacteria begin 
    innerjoin(tidy_spot, on = :B_ID)
    transform(
        [:M_X, :B_X] => ((m, b) -> m - b) => :RELATIVE_SPOT_X, 
        [:M_Y, :B_Y] => ((m, b) -> m - b) => :RELATIVE_SPOT_Y)
    groupby([:T_ID, :FRAME])
    combine(
        :N_MAXIMA => ((x) -> allequal(x)[2]) => :N_MAXIMA,
        :FRAME_START => ((x) -> allequal(x)[2]) => :FRAME_START,
        :RELATIVE_SPOT_X => ((x) -> vec([x])) => :RELATIVE_SPOT_X,
        :RELATIVE_SPOT_Y => ((x) -> vec([x])) => :RELATIVE_SPOT_Y)
    groupby([:T_ID])
    transform(
        :RELATIVE_SPOT_X => ((x) -> lag(x, default = [0])) => :PREVIOUS_SPOT_X,
        :RELATIVE_SPOT_Y => ((x) -> lag(x, default = [0])) => :PREVIOUS_SPOT_Y,
        :RELATIVE_SPOT_X => ((x) -> lag(x, 2, default = [0])) => :PREVIOUS2_SPOT_X,
        :RELATIVE_SPOT_Y => ((x) -> lag(x, 2, default = [0])) => :PREVIOUS2_SPOT_Y,
        :RELATIVE_SPOT_X => ((x) -> lead(x, default = [0])) => :NEXT_SPOT_X,
        :RELATIVE_SPOT_Y => ((x) -> lead(x, default = [0])) => :NEXT_SPOT_Y)
end



df[!, :RELATIVE_SPOT_X] = [(length(x₋₁) > length(x) < length(x₊₁)) ? mean(cat(x₋₁, x₋₂, dims = 2), dims = 2) : x for (x, x₋₁, x₋₂, x₊₁) in zip(df[!, :RELATIVE_SPOT_X], df[!, :PREVIOUS_SPOT_X], df[!, :PREVIOUS2_SPOT_X], df[!, :NEXT_SPOT_X])]

df[!, :ΔX] = [ismissing(x₋₁) ? missing : vec(x' .- x₋₁) for (x, x₋₁) in zip(df[!, :RELATIVE_SPOT_X], df[!, :PREVIOUS_SPOT_X])]
df[!, :ΔY] = [ismissing(y₋₁) ? missing : vec(y' .- y₋₁) for (y, y₋₁) in zip(df[!, :RELATIVE_SPOT_Y], df[!, :PREVIOUS_SPOT_Y])]
##
df[!, :ΔR] = [radialcoordinate.(Δx', Δy)  for (Δx, Δy) in zip(df[!, :ΔX], df[!, :ΔY])]
##
df[!, :CLOSEST_FOCIS] = [ismissing(Δr) ? missing : argmin(Δr, dims = 1) for Δr in df[!, :ΔR]]
##
df[!, :ΔR] = [ismissing(Δr) ? missing : Δr[i] for (Δr, i) in zip(df[!, :ΔR], df[!, :CLOSEST_FOCIS])]
##
df[!, :ΔX] = [ismissing(Δx) ? missing : Δx[:,:]'[i] for (Δx, i) in zip(df[!, :ΔX], df[!, :CLOSEST_FOCIS])]
##
df[!, :ΔY] = [ismissing(Δy) ? missing : Δy[:,:][i] for (Δy, i) in zip(df[!, :ΔY], df[!, :CLOSEST_FOCIS])]
##

df[!, [:ΔX]] = [Δx[sortperm(Δr)] for (Δx, Δy, Δr) in zip(df[!, :ΔX], df[!, :ΔY], df[!, :ΔR])]

df[!, :θ]  = [angularcoordinate.(Δx, Δy) for (Δx, Δy) in zip(df[!, :ΔX], df[!, :ΔY])]
df[!, :CLOSEST_FOCIS] = [ismissing(Δr) ? missing : argmin(Δr, dims = 1) for Δr in df[!, :ΔR]]
df[!, :ΔR] = [ismissing(Δr) ? missing : Δr[i] for (Δr, i) in zip(df[!, :ΔR], df[!, :CLOSEST_FOCIS])]
df[!, :ΔX] = [ismissing(Δx) ? missing : Δx[i] for (Δx, i) in zip(df[!, :ΔX], df[!, :CLOSEST_FOCIS])]
df[!, :ΔY] = [ismissing(Δy) ? missing : Δy[i] for (Δy, i) in zip(df[!, :ΔY], df[!, :CLOSEST_FOCIS])]
df[!, :θ]  = [ismissing(θ)  ? missing : θ[i]  for (θ, i)  in zip(df[!, :θ],  df[!, :CLOSEST_FOCIS])]

df = transform(df, :θ => lag => :PREVIOUS_θ)

df[!, :OPPOSITE] = [ismissing(θ) || ismissing(θ₋₁) ? missing : isoppositedirection.(θ', θ₋₁, 45) for (θ, θ₋₁) in zip(df[!, :θ], df[!, :PREVIOUS_θ])]

df[!, :ΔR] = [ismissing(opp) ? Δr : opp'[i].*Δr for (opp, i, Δr) in zip(df[!, :OPPOSITE], df[!, :CLOSEST_FOCIS], df[!, :ΔR])]
df

##
test = @chain df begin
    filter(row -> row.T_ID == "tb11", _)
    select([:T_ID, :FRAME, :N_MAXIMA, :FRAME_START, :ΔX, :ΔY, :ΔR, :θ])
#    transform(:ΔR => cumsum => :MOVEMENT)
#    flatten([:ΔX, :ΔY, :ΔR, :θ])
end
##

#test[!, :ZSCORE] = test[!, :ZSCORE] / 10
set_default_plot_size(50cm, 30cm)
@pipe test |>
  #  filter(row -> row.ZSCORE >= 1, _) |>
    plot(_,
        x=:FRAME,
        y=:MOVEMENT,
  #      size=:RELATIVE_INTENSITY,
        Geom.point,
        Guide.xticks(ticks=1:40)) |> 
        PNG("foo.png", 50cm, 30cm)


function cross_combinations_difference(row_to_cross, col_to_cross)
    m = Array{Float64}(undef, length(row_to_cross), length(col_to_cross))
    for (col_i, col_v) in enumerate(col_to_cross)
        for (row_i, row_v) in enumerate(row_to_cross)
            m[row_i, col_i] = row_v - col_v
        end
    end
    return m
end

minimum(filter(!isnan, df.ANGLE))

@pipe df |> 
    sort(_,
    [:FRAME, :ANGLE]) |>
    groupby(_, 
    :T_ID)