using ImageNoise
using BioformatsLoader
using ImageView
using JavaCall
using PyCall
using FileIO
using ImageAxes
using SimpleTraits
using Wavelets
using BM3D
using ImagePhaseCongruency
using PerceptualColourMaps
using Colors
using Plots
using Images
using TiffImages
@traitimpl TimeAxis{Axis{:t}}


BioformatsLoader_path = "/home/thomas/BioinformaticsTools/bioformats_package.jar"
JavaCall.addClassPath(BioformatsLoader_path)
JavaCall.addOpts("-ea") # Enable assertions
JavaCall.addOpts("-Xmx28672M")

try
    JavaCall.init()
catch
    println("Failed to initialize Bioformat")
end

img_path = "/run/media/thomas/DATA/MICROSCOPY/210504-GL1205-2h-PI.nd2"
img = bf_import(img_path)[1]
img = permutedims(img,[2, 5, 3, 4, 1])
img = AxisArray(img, :r, :c, :x, :y, :t)

img

size(img[c = 1, r = 1])
ROIs = size(img)[2]


for time in 1:timepoints
    save("~/img_$time.",  img[t = time, r = 1])
end

temp = data(data(img[c = 1, r = 1]))

axisnames(img)
data(img[t = 5, r = 1])
arraydata()
TiffImages.DenseTaggedImage(temp)

save("img_5.tiff",  data(data(img[r = 1, c = 1])))
size(img)
display(data(img)["Pixels"])
imshow(frame)
save("~/img.ome.tiff",  img[t = 1, r = 1, c = 2])
save("~/img_denoised.tiff",  floor.(UInt16,denoised))
save("~/img_denoised_refilter.tiff",  floor.(UInt16,denoised_refilter))
save("~/img_denoised_nlm.tiff",  floor.(UInt16,denoised_refilter))

convert(Array{UInt16} , denoised)

imshow(floor.(UInt16,denoised))


permutedims

frame = img[t = 1, c = 2, r = 1]

(pc, or, ft, T) =
         phasecongmono(ppdrc(Float64.(frame)/1e5, 50); nscale=4, minwavelength=3, mult=2,
                        sigmaonf=0.55, k=3, cutoff=0.5, g=10,
                        deviationgain=1.5, noisemethod=-1)
(phaseSym, symmetryEnergy, T) = phasesymmono(ppdrc(Float64.(frame)/1e5, 50); nscale=5, polarity=1)
imshow(phaseSym)

heatmap(pc, fill  = ColorGradient(cmap("L16")))
plot!(size=(2048,2048))


heatmap(phaseSym, fill  = ColorGradient(cmap("L16")))
plot!(size=(2048,2048))

heatmap(frame, fill  = ColorGradient(cmap("L16")))
plot!(size=(2048,2048))

scale = 1e4
imshow(ppdrc(Float64.(frame)/1e5, 25))
imshow(frame)

imshow(histeq(frame, 100000))

adjusted = adjust_histogram(Float64.(frame)/1e5, Equalization(nbins = 1e7))
imshow(adjusted)

Float64.(frame)
img_t = Float64.(testimage("m51"))
adjusted_t = adjust_histogram(Float64.(img_t), Equalization(nbins = 100))
imshow(adjusted_t)
imshow(img_t)


median(Float64.(frame)/1e5)
median(img_t)

imshow(pc)
(P, S, p, s) = perfft2(Float64.(frame)) 
imshow(s)