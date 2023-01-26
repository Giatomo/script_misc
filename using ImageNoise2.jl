using ImageNoise
using BioformatsLoader
using ImageView
using JavaCall
using PyCall
using FileIO
using SimpleTraits
using Wavelets
using BM3D
using Images
using ImageFiltering
using PaddedViews
import ImageCore.channelview
import Base.sign
using TestImages
using Colors
@traitimpl TimeAxis{Axis{:t}}


img = testimage("blobs") # fullname
img = Gray.(img)

BioformatsLoader_path = "/home/thomas/BioinformaticsTools/bioformats_package.jar"
JavaCall.addClassPath(BioformatsLoader_path)
JavaCall.addOpts("-ea") # Enable assertions
JavaCall.addOpts("-Xmx8192M")

try
    JavaCall.init()
catch
    println("Failed to initialize Bioformat")
end

img_path = "/run/media/thomas/DATA/MICROSCOPY/20210209_ZapTGFP_TC_with_TL_between/3h001.nd2"
image = bf_import(img_path)[1]
image = AxisArray(image, :t, :r, :x, :y, :c)

img = image[t = 1, r = 1, c = 2]

function channelview(A::Float64)::Float64
    return A
end


function chain(from, to::Matrix{CartesianIndex{2}}, Ɛₓ, Ɛᵧ, gₓ, gᵧ, half_search_window)

    @show from 
    Δₓ =  Ɛₓ[to] .- Ɛₓ[from]
    Δᵧ =  Ɛᵧ[to] .- Ɛᵧ[from] 
    Δᵣ = .√(Δₓ.^2+Δᵧ.^2)

    score = sign.(gᵧ[from] .* Δₓ .- gₓ[from] .* Δᵧ) ./ Δᵣ
    invalid_pos = ((gᵧ[from] .* Δₓ .- gₓ[from] .* Δᵧ) .* (gᵧ[to] .* Δₓ .- gₓ[to] .* Δᵧ )) .<= 0

    if invalid_pos == false
        return score[1]
    end
    if invalid_pos == true
        return [0]
    end

    limits = Tuple(size(Ɛₓ)) .- half_search_window
    
    @show (Tuple.(to))
    invalid_pos = invalid_pos .| map(x -> any(x .> limits), Tuple.(to)) 
    score[invalid_pos] .= 0
    return score

end

function chain(from, to::CartesianIndex{2}, Ɛₓ, Ɛᵧ, gₓ, gᵧ, half_search_window)

    @show from 
    Δₓ =  Ɛₓ[to] .- Ɛₓ[from]
    Δᵧ =  Ɛᵧ[to] .- Ɛᵧ[from] 
    Δᵣ = .√(Δₓ.^2 + Δᵧ.^2)

    score = sign.(gᵧ[from] .* Δₓ .- gₓ[from] .* Δᵧ) ./ Δᵣ
    invalid_pos = ((gᵧ[from] .* Δₓ .- gₓ[from] .* Δᵧ) .* (gᵧ[to] .* Δₓ .- gₓ[to] .* Δᵧ )) .<= 0

    if invalid_pos == false
        return score[1]
    end
    if invalid_pos == true
        return [0]
    end

    limits = Tuple(size(Ɛₓ)) .- half_search_window
    invalid_pos = invalid_pos | any(Tuple(to) > limits) 
    score[invalid_pos] .= 0
    return score

end
(prev, next) = chain_edge_points(Ɛₓ, Ɛᵧ, gₓ, gᵧ)

function image_gradient(img, scale = 3)
    img_s = imfilter(img, Kernel.gaussian(scale))

    gₓ = img_s[:, 1:end-2]  .- img_s[:, 3:end]
    gᵧ = img_s[1:end-2, :] .-  img_s[3:end, :]

    gₓ = PaddedView(0, gₓ, (size(gₓ)[1],size(gₓ)[2]+2), (1,2))
    gᵧ = PaddedView(0, gᵧ, (size(gᵧ)[1]+2,size(gᵧ)[2]), (2,1))

    g  = .√(gₓ.^2 .+ gᵧ.^2)
    return (g, gₓ, gᵧ)
end

(g, gₓ, gᵧ) = image_gradient(img)

function compute_edge_points(g, gₓ, gᵧ)
    θₓ = g[:, 1:end-2] .> g[:, 2:end-1]  .>= g[:, 3:end]
    θₓ = PaddedView(0, θₓ, (size(θₓ)[1],size(θₓ)[2]+2), (1,2))
    θₓ =  θₓ .& (gₓ .>= gᵧ)

    θᵧ = g[1:end-2, :] .> g[2:end-1, :]  .>= g[3:end, :]
    θᵧ = PaddedView(0, θᵧ, (size(θᵧ)[1]+2,size(θᵧ)[2]), (2,1))
    θᵧ = θᵧ .& (gₓ .<= gᵧ)

    Aₓ  = g[:, 1:end-2]
    Aₓ  = PaddedView(0, Aₓ, (size(Aₓ)[1],size(Aₓ)[2]+2), (1,2))
    Aᵧ  = g[1:end-2, :]
    Aᵧ  = PaddedView(0, Aᵧ, (size(Aᵧ)[1]+2,size(Aᵧ)[2]), (2,1))
    Aₓᵧ = g[1:end-2, 1:end-2]
    Aₓᵧ = PaddedView(0, Aₓᵧ, (size(Aₓᵧ)[1]+2,size(Aₓᵧ)[2]+2), (2,2))

    Cₓ  = g[:, 3:end]
    Cₓ  = PaddedView(0, Cₓ, (size(Cₓ)[1],size(Cₓ)[2]+2))
    Cᵧ  = g[3:end, :]
    Cᵧ  = PaddedView(0, Cᵧ, (size(Cᵧ)[1]+2,size(Cᵧ)[2]))
    Cₓy = g[3:end, 3:end]
    Cₓy = PaddedView(0, Cₓy, (size(Cₓy)[1]+2,size(Cₓy)[2]+2))

    λₓ  = (Aₓ .- Cₓ) ./ 2(Aₓ .- 2g .+ Cₓ)
    λᵧ  = (Aᵧ .- Cᵧ) ./ 2(Aᵧ .- 2g .+ Cᵧ)
    λₓᵧ = (Aₓᵧ .- Cₓy) ./ 2(Aₓᵧ .- 2g .+ Cₓy)

    x_ind = (1:size(g)[2])' .* ones(size(g))
    y_ind = (1:size(g)[1])  .* ones(size(g))

    λ =  (θₓ .& θᵧ) .* λₓᵧ .+ (θₓ .> θᵧ) .* λₓ .+ (θᵧ .> θₓ) .* λᵧ
    Ɛₓ = (θₓ .| θᵧ) .* x_ind + λ .* θₓ
    Ɛₓ[.~(θₓ .| θᵧ)] .= -1
    Ɛᵧ = (θₓ .| θᵧ) .* y_ind + λ .* θᵧ
    Ɛᵧ[.~(θₓ .| θᵧ)] .= -1
    
    #Ɛ = [ϵ for ϵ in collect(zip(Ɛₓ,Ɛᵧ)) if ϵ != (0,0)]
    return (Ɛₓ, Ɛᵧ)
end

(Ɛₓ, Ɛᵧ) = compute_edge_points(g, gₓ, gᵧ)

origin = CartesianIndex(0, 0)
out = CartesianIndex(-1, -1)

function chain_edge_points(Ɛₓ, Ɛᵧ, gₓ, gᵧ, search_window = (5,5))
    half_search_window = floor.(Int, search_window./2)
    prev = next = fill(out, size(Ɛₓ))
    for from in CartesianIndices(Ɛₓ)[1:end-half_search_window[1], 1:end-half_search_window[2]]
        if (Ɛₓ[from] <= 0) | (Ɛᵧ[from] <= 0)
            continue
        end
        try
            to = from - CartesianIndex(half_search_window[1], half_search_window[2]): from + CartesianIndex(half_search_window[1], half_search_window[2])
            score = chain(from, to, Ɛₓ, Ɛᵧ, gₓ, gᵧ, half_search_window)
            global (fwd_score, fwd_pos) = findmax(score)
            fwd_pos += from - CartesianIndex(half_search_window[1], half_search_window[2])
            global (bck_score, bck_pos) = findmin(score)
            bck_pos += from - CartesianIndex(half_search_window[1], half_search_window[2])
            

        catch
            continue
        end

        if (fwd_pos >= origin) && 
           (next[from] != fwd_pos) && 
           (((alt = prev[fwd_pos]) < origin) || all(chain(alt, fwd_pos,  Ɛₓ, Ɛᵧ, gₓ, gᵧ, half_search_window) .< fwd_score))

           if next[from] >= origin 
                prev[next[from]] = out
            end
            next[from] = fwd_pos
            if alt >= origin
                next[alt] = out
            end
            prev[fwd_pos] = from
        end

        if (bck_pos >= origin) && 
           (prev[from] != bck_pos) && 
           (((alt = next[bck_pos])) < origin) || all(chain(alt, bck_pos,  Ɛₓ, Ɛᵧ, gₓ, gᵧ, half_search_window) .> bck_score)

            if alt >= origin
                prev[alt] = out
            end
            next[bck_pos] = from
            if prev[from] >= origin
                next[prev[from]] = out
            end
            prev[from] = bck_pos
        end
    end
    
    return (prev, next)
end

(prev, next) = chain_edge_points(Ɛₓ, Ɛᵧ, gₓ, gᵧ)


minimum(next)

CartesianIndex(154,189) -CartesianIndex(2,2) : CartesianIndex(154,189) + CartesianIndex(2,2)



σ = Threshold.noisest(img[t = 1, r = 1, c = 2])
denoised = floor.(UInt16, bm3d(convert(Array{Float64},img[t = 1, r = 1, c = 2]), σ))

imgedg = canny(denoised, (Percentile(1), Percentile(99)))

imshow(imgedg)

σ = Threshold.noisest(img[t = 1, r = 1, c = 2])
σ = Threshold.noisest(denoised)





(to, limit) = (CartesianIndex{2}[CartesianIndex(233, 252) CartesianIndex(233, 1000) CartesianIndex(1000, 254) CartesianIndex(233, 255) CartesianIndex(233, 256); CartesianIndex(234, 252) CartesianIndex(234, 253) CartesianIndex(234, 254) CartesianIndex(234, 255) CartesianIndex(234, 256); CartesianIndex(235, 252) CartesianIndex(235, 253) CartesianIndex(235, 254) CartesianIndex(235, 255) CartesianIndex(235, 256); CartesianIndex(236, 252) CartesianIndex(236, 253) CartesianIndex(236, 254) CartesianIndex(236, 255) CartesianIndex(236, 256); CartesianIndex(237, 252) CartesianIndex(237, 253) CartesianIndex(237, 254) CartesianIndex(237, 255) CartesianIndex(237, 256)], CartesianIndex(254, 256))
for coord in Tuple.(to)
    @show all(coord .<= Tuple(limit))
end

map(x -> all(x .> Tuple(limit)), Tuple.(to)) 
map(x -> all(x .> Tuple(limit)), Tuple.(to[1])) 
to[1]

Tuple(to)
length(to)
length(limit)
typeof(to)
typeof(limit)

(233, 252) < (254, 256)
(233, 258) < (254, 256)
(255, 258) < (254, 256)
all((255, 252) .< (254, 256))




Tuple(CartesianIndex(0, 2)) .<= [Tuple(CartesianIndex(1, 1))] # T

CartesianIndex(1, 0) <= CartesianIndex(1, 1)  # T

CartesianIndex(0, 1) <= CartesianIndex(1, 1)  # T

CartesianIndex(1, 1) <= CartesianIndex(1, 1)  # T
CartesianIndex(2, 1) <= CartesianIndex(1, 1)  # F
CartesianIndex(1, 2) <= CartesianIndex(1, 1)  # F
CartesianIndex(1, 2) <= CartesianIndex(2, 1)  # F


CartesianIndex(255, 255) <= CartesianIndex(254, 256)

tmp = Gray{Float64}[Gray{Float64}(-0.029369527889321054) Gray{Float64}(-0.0465157433727421) Gray{Float64}(-0.08626558697677589) Gray{Float64}(-0.05828726043015044) Gray{Float64}(-0.05984265442736327); Gray{Float64}(-0.008488231619196342) Gray{Float64}(-0.016505012396733108) Gray{Float64}(-0.05176164600091732) Gray{Float64}(0.11030392815893038) Gray{Float64}(-4.985453074856153); Gray{Float64}(0.010820758030480982) Gray{Float64}(0.007411307570152676) Gray{Float64}(0.0) Gray{Float64}(0.029515908161280427) Gray{Float64}(0.04975953980623842); Gray{Float64}(0.028364716525147778) Gray{Float64}(0.026020321635306688) Gray{Float64}(0.023255380184090523) Gray{Float64}(0.0416839094783496) Gray{Float64}(-4.1961181036444595); Gray{Float64}(0.04554603584846081) Gray{Float64}(0.04327486883538903) Gray{Float64}(0.04122133456305061) Gray{Float64}(0.060619582242914416) Gray{Float64}(-4.1961181036444595)]
print(tmp)
sign(tmp)

sign(Gray{Float64}(-0.029369527889321054))
sign.(tmp)

sign(Gray{Float64}(-0.029369527889321054))
print(tmp)[1]
typeof(tmp[1])
function sign(A::Gray{Float64})::Float64
    return sign(Float64(A))
end

sign.(Gray{Float64}(-0.02452245386462))