from pathlib import Path
import numpy as np
import tifffile as tf
from pims import ND2Reader_SDK
from skimage.io import imread, imsave
from skimage.restoration import denoise_nl_means, estimate_sigma
from bm3d import bm3d, BM3DProfile
from tqdm import tqdm

# Directory select and create
WORK_PATH = Path("/run/media/thomas/DATA/MICROSCOPY")
DATA_PATH = WORK_PATH/"20210209_ZapTGFP_TC_with_TL_between"

FILE = DATA_PATH/"3h002.nd2"
PHASE_CHANNEL = 2
BDELLO_FLUO_CHANNEL = 1 


FramesSequenceND.__init__()
with ND2Reader_SDK(FILE) as img:

    img.bundle_axes = "txyc"
    img.iter_axes = "m"

# ---------------------------- Image metadata ---------------------------- #

    width      = img.sizes["x"]
    heigth     = img.sizes["y"]
    n_roi      = img.sizes["m"]
    channels   = img.sizes["c"]
    timepoints = img.sizes["t"]

    px_size    = img.metadata["calibration_um"]

    print(np.shape(img[1]))
# ---------------------------------------------------------------------------- #


    for roi in tqdm(range(n_roi)):

        for t in tqdm(range(timepoints)):

            frame = img[roi][t,:,:,1]
            print()
            frame = np.ascontiguousarray(frame)

            σ = np.mean(estimate_sigma(frame, multichannel = False))

            img[roi][t,:,:,1] = bm3d(frame, σ)


tf.imwrite(f'denoised{roi}.tif', temp_t,
        resolution = (1./px_size, 1./px_size),
        metadata = {'unit': 'um', 'axes': 'XYT'})




plt.subplot(1,2,1)
plt.imshow(np.asarray(img_denoised),  plt.get_cmap("Greys_r"))
plt.subplot(1,2,2)
plt.imshow(np.asarray(img), plt.get_cmap("Greys_r"))
plt.show()




tf.imwrite('denoised.tif', temp_roi,
        resolution = (1./px_size, 1./px_size),
        metadata = {'unit': 'um', 'axes': 'XYTM'})



with ND2Reader_SDK(FILE) as frames:

    frames.bundle_axes = "tmxyc"
    print(frames.shape)


