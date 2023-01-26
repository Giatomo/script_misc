import random as rd
import shutil
from itertools import chain
from pathlib import Path

import cv2
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PIL
from PIL import Image
from PIL import Image as im
from PIL import ImageOps
from pims import ND2Reader_SDK
from skimage.exposure import rescale_intensity
from skimage.filters import threshold_local, threshold_yen
from skimage.io import imread, imsave
from tensorflow import keras
from tensorflow.keras import layers
from tqdm import tqdm




# Directory select and create
WORK_PATH = Path(Path.home()/"CNN"/"Bdello_UNet")
DATA_PATH = WORK_PATH/"images"
SEGMENTATION_PATH = WORK_PATH/"segmentation"
BDELLO_FLUO_PATH = WORK_PATH/"bd_fluo"
BDELLO_FLUO_YEN_PATH = WORK_PATH/"bd_fluo_yen"

TL_WORK_PATH = Path(Path.home()/"CNN"/"Bdello_UNet_correction")
TL_DATA_PATH = TL_WORK_PATH/"images"
TL_SEGMENTATION_PATH = TL_WORK_PATH/"segmentation"
TL_BDELLO_FLUO_PATH = TL_WORK_PATH/"bd_fluo"
TL_BDELLO_FLUO_YEN_PATH = TL_WORK_PATH/"bd_fluo_yen"

INPUT_TL_FILES = list(Path("/mnt/Renske/Renske - TdTomato growth").glob('**/*.tif'))

folder_to_create = [DATA_PATH, SEGMENTATION_PATH, BDELLO_FLUO_PATH, BDELLO_FLUO_YEN_PATH, TL_DATA_PATH, TL_SEGMENTATION_PATH, TL_BDELLO_FLUO_PATH, TL_BDELLO_FLUO_YEN_PATH]
for folder in folder_to_create:
    folder.mkdir(parents = True, exist_ok = True)

IMG_SIZE     = (128,128,1)
PADDING      = int(IMG_SIZE[0]/2)

PHASE_CHANNEL = 2
BDELLO_FLUO_CHANNEL = 1 

# Load ND2 file with the Bdelloplast 
INPUT_FILES = Path("//home/thomas/Bureau/2021-01-12_TL_training_data_CNN").glob('**/*.nd2')
img_i = 0

# Load .nd2 file and threshold the fluorescence to use as the "correct" segmentation for phase
for nd2_file in INPUT_FILES:

    csv_file = nd2_file.with_suffix(".csv")

    with ND2Reader_SDK(nd2_file) as frames:

        N_LOCATIONS = frames.sizes["m"]
        PIXEL_RATIO = frames.metadata["calibration_um"]

        frames.iter_axes = "m"
        frames.bundle_axes = "xyc"


        # To access one image: frames[time][position, X, Y, channel]
        cell_positions = pd.read_csv(csv_file, sep=";")
        cell_positions[["X", "Y"]] = (cell_positions[["LOCATION.x", "LOCATION.y"]]/PIXEL_RATIO).astype(int)
        cell_positions[["FILENAME", "IMG_N"]] = cell_positions["IMAGE.meta"].str.split(r"c:\d+\/\d+ - (\w+\.\w+) \(series (\d+)\)", expand=True)[[1,2]]
        cell_positions = cell_positions[["X", "Y", "FILENAME", "IMG_N"]] 
        cell_positions["IMG_N"] = cell_positions["IMG_N"].astype(int)

        N_CELLS = len(cell_positions)
        # Generating training dataset
        for loc_i in range(N_LOCATIONS):

            img = np.flipud(np.rot90(frames[loc_i]))

            retval, binary_mask = cv2.threshold(img[:, :, BDELLO_FLUO_CHANNEL], 0, 1, cv2.THRESH_OTSU)

            temp = img[:, :, BDELLO_FLUO_CHANNEL]
            temp[binary_mask == 0] = 0
            block_size = 11
            local_thresh = threshold_local(temp, block_size)
            bin_local = np.asarray(temp > local_thresh, dtype = "uint16")
            
            for x, y in zip(cell_positions[cell_positions["IMG_N"] == loc_i+1]["X"], cell_positions[cell_positions["IMG_N"] == loc_i+1]["Y"]):

                try :
                    # Cropping around bdelloplast in phase -> INPUT
                    data = im.fromarray(img[y-PADDING:y+PADDING, x-PADDING:x+PADDING, PHASE_CHANNEL])
                    data.save(DATA_PATH/f"img_{img_i:04}.png")

                    # Cropping around bdelloplast in fluorescence + making binary mask -> SEGMENTATION 
                    segmentation = im.fromarray(bin_local[y-PADDING:y+PADDING, x-PADDING:x+PADDING])
                    segmentation.save(SEGMENTATION_PATH/f"img_{img_i:04}.png")

                    bd_fluo = im.fromarray(img[y-PADDING:y+PADDING, x-PADDING:x+PADDING, BDELLO_FLUO_CHANNEL])
                    bd_fluo.save(BDELLO_FLUO_PATH/f"img_{img_i:04}.png")

                    img_i += 1

                except(SystemError):
                    continue



    csv_tl_file = Path("/home/thomas/bdello_pos.csv")

    
    from tifffile import imread, imwrite, TiffFile, TiffWriter

    cell_positions = pd.read_csv(csv_tl_file, sep=";")
    cell_positions[["X", "Y"]] = (cell_positions[["LOCATION.x", "LOCATION.y"]]/0.07).astype(int)
    cell_positions = cell_positions[["X", "Y"]] 

    img_t = 0
    for last_phase, last_fluo, phase, fluo in zip(tqdm(TiffFile(INPUT_TL_FILES[0]).pages[1:]), TiffFile(INPUT_TL_FILES[1]).pages[1:], list(TiffFile(INPUT_TL_FILES[0]).pages[:-1]), TiffFile(INPUT_TL_FILES[1]).pages[:-1]):
        img_t += 1
        
        
        phase = phase.asarray()
        fluo = fluo.asarray()
        last_phase = last_phase.asarray()
        last_fluo = last_fluo.asarray()


        # phase = np.flipud(np.rot90(phase))
        # fluo = np.flipud(np.rot90(fluo))
        
        retval, binary_mask = cv2.threshold(fluo, 0, 1, cv2.THRESH_OTSU)
        retval, last_binary_mask = cv2.threshold(last_fluo, 0, 1, cv2.THRESH_OTSU)

        temp = fluo
        temp[binary_mask == 0] = 0
        block_size = 11
        local_thresh = threshold_local(temp, block_size)
        bin_local = np.asarray(temp > local_thresh, dtype = "uint16")



        temp = last_fluo
        temp[binary_mask == 0] = 0
        block_size = 11
        local_thresh = threshold_local(temp, block_size)
        last_bin_local = np.asarray(temp > local_thresh, dtype = "uint16")




        img_i = 0
        for x, y in zip(cell_positions["X"], cell_positions["Y"]):
            print(img_i)
            if x < PADDING or y < PADDING or x + PADDING > 1608 or y + PADDING > 1608:
                print("skiping...")
                continue
            img_i += 1

            try :
                # Cropping around bdelloplast in phase -> INPUT
                data = np.stack([last_bin_local[y-PADDING:y+PADDING, x-PADDING:x+PADDING], phase[y-PADDING:y+PADDING, x-PADDING:x+PADDING]], axis=0)
                imwrite(TL_DATA_PATH/f"t_{img_t:04}_img_{img_i:04}.tiff", data)


                # Cropping around bdelloplast in fluorescence + making binary mask -> SEGMENTATION 
                segmentation = bin_local[y-PADDING:y+PADDING, x-PADDING:x+PADDING]
                imwrite(TL_SEGMENTATION_PATH/f"t_{img_t:04}_img_{img_i:04}.tiff", segmentation)

            except(SystemError):
                continue

log = [] 
def remove_data(data_path, segm_path):
    os.remove(data_path)
    os.remove(segm_path)
    print(data_path)
    log.append(data_path)
    plt.close("all")
    return 


from matplotlib.widgets import Button
import os

ALL_DATA = sorted(list(TL_DATA_PATH.glob('**/*.tiff')))
ALL_SEGM = sorted(list(TL_SEGMENTATION_PATH.glob('**/*.tiff')))

for data, segm in zip(ALL_DATA, ALL_SEGM):
    img = imread(data)
    last_binary = img[0,:,:]
    phase = img[1,:,:]
    binary = imread(segm)

    plt.subplot(1,3,1)

    plt.imshow(phase, cmap= plt.get_cmap("Greys_r"))
    plt.axis('off')
    plt.title('Phase')
    plt.subplot(1,3,2)
    plt.imshow(last_binary, plt.get_cmap("binary_r"))
    plt.axis('off')
    plt.title('Last segmentation')
    plt.subplot(1,3,3)
    plt.imshow(binary, plt.get_cmap("binary_r"))
    plt.axis('off')
    plt.title('Segmentation')
    axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
    axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
    bnext = Button(axnext, 'Delete')
    bnext.on_clicked(lambda x : remove_data(data, segm))
    bprev = Button(axprev, 'Keep')
    bprev.on_clicked(lambda x: plt.close("all"))
    plt.show()

Counter((re.match(r'.*img_(\d{4}).*', str(i))[1] for i in log))

print([str(i) for i in log])
print(re.match(r'.*img_(\d{4}).*', str(log[0]))[1])

log = []


ALL_DATA = sorted(list(TL_DATA_PATH.glob('**/*.tiff')))
ALL_SEGM = sorted(list(TL_SEGMENTATION_PATH.glob('**/*.tiff')))

common = list(set([p.name for p in ALL_SEGM]) & set([p.name for p in ALL_DATA]))

rows = [[data, re.match(r'.*img_(\d{4}).*', str(data))[1], int(re.match(r'.*t_(\d{4}).*', str(data))[1]), np.median(imread(data)[1,:,:]), np.std(imread(data)[1,:,:])] for data in ALL_DATA]

df = pd.DataFrame(rows, columns=["path", "img", "t", "value", "value_sd"])
df = df.groupby('img').apply(lambda grp: grp.assign(mean=grp['value'].mean()))
df = df.groupby('img').apply(lambda grp: grp.assign(sd=grp['value'].std()))

filtered = df[~df["value"].between(df["mean"] - 2*df["sd"], df["mean"] + 2*df["sd"])]
filtered = filtered[filtered["t"] >= 10]
filtered.groupby("img").min()
sns.relplot(data=filtered, x='t', y = "value", hue='img', kind = "line")
plt.show()

log = []
for p in filtered.groupby("img").min()["path"]:
    plt.imshow(imread(p)[1,:,:], cmap= plt.get_cmap("Greys_r"))
    axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
    axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
    bnext = Button(axnext, 'Delete')
    bnext.on_clicked(lambda x : [log.append(p), plt.close("all")])
    axprev = Button(axprev, 'Keep')
    axprev.on_clicked(lambda x : [ plt.close("all")])
    plt.show()



log = []
ALL_DATA = sorted(list(TL_DATA_PATH.glob('**/*.tiff')))
ALL_SEGM = sorted(list(TL_SEGMENTATION_PATH.glob('**/*.tiff')))
for data, segm in zip(ALL_DATA, ALL_SEGM):
    img = imread(data)
    last_binary = img[0,:,:]
    phase = img[1,:,:]
    binary = imread(segm)

    if ( int(re.match(r'.*t_(\d{4}).*', str(data))[1]) > 20):
        plt.subplot(1,3,1)
        plt.imshow(phase, cmap= plt.get_cmap("Greys_r"))
        plt.axis('off')
        plt.title('Phase')
        plt.subplot(1,3,2)
        plt.imshow(last_binary, plt.get_cmap("binary_r"))
        plt.axis('off')
        plt.title('Last segmentation')
        plt.subplot(1,3,3)
        plt.imshow(binary, plt.get_cmap("binary_r"))
        plt.axis('off')
        plt.title('Segmentation')
        axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
        axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
        bnext = Button(axnext, 'Delete')
        bnext.on_clicked(lambda x : remove_data(data, segm))
        bprev = Button(axprev, 'Keep')
        bprev.on_clicked(lambda x: plt.close("all"))
        plt.show()
        print(sorted(log))

print(Counter((re.match(r'.*img_(\d{4}).*', str(i))[1] for i in log)))
print(log)
# ALL_FLUO = sorted(list(BDELLO_FLUO_PATH.glob('**/*.png')))

for data, segm in zip(ALL_DATA[:10], ALL_SEGM[:10]):
    print(data, "|", segm)

# Split our img paths into a training and a validation set

seed = rd.randint(0, 1e100)
val_samples = int(len(ALL_DATA)*0.7)
rd.Random(seed).shuffle(ALL_DATA)
rd.Random(seed).shuffle(ALL_SEGM)
# rd.Random(seed).shuffle(ALL_FLUO)
train_input_img_paths = ALL_DATA[:-val_samples]
train_target_img_paths = ALL_SEGM[:-val_samples]
val_input_img_paths = ALL_DATA[-val_samples:]
val_target_img_paths = ALL_SEGM[-val_samples:]
# val_fluo_img_paths = ALL_FLUO[-val_samples:]


from tensorflow.keras.preprocessing.image import load_img

class TrainingData(keras.utils.Sequence):
    """Helper to iterate over the data (as Numpy arrays)."""

    def __init__(self, batch_size, img_size, input_img_paths, target_img_paths):
        self.batch_size = batch_size
        self.img_size = img_size
        self.input_img_paths = input_img_paths
        self.target_img_paths = target_img_paths


    def __len__(self):
        return len(self.target_img_paths) // self.batch_size

    def __getitem__(self, idx):
        """Returns tuple (input, target) correspond to batch #idx."""
        i = idx * self.batch_size
        batch_input_img_paths = self.input_img_paths[i : i + self.batch_size]
        batch_target_img_paths = self.target_img_paths[i : i + self.batch_size]
        x = np.zeros((self.batch_size,) + self.img_size + (2,), dtype="float32")
        for j, path in enumerate(batch_input_img_paths):
            # img = load_img(path, target_size=self.img_size,  color_mode=self.color_mode)
            img = imread(path)
            img = np.moveaxis(img, 0, 2)
            x[j] = img
        y = np.zeros((self.batch_size,) + self.img_size + (1,), dtype="uint8")
        for j, path in enumerate(batch_target_img_paths):
            # img = load_img(path, target_size=self.img_size, color_mode=self.color_mode)
            img = imread(path)
            img = np.expand_dims(img, 2)
            # print("y :", img.shape)
            y[j] = img
            # Ground truth labels are 1, 2, 3. Subtract one to make them 0, 1, 2:
        return x, y


def get_model(img_size, num_classes):
    inputs = keras.Input(shape=img_size + (2,))


    ### [First half of the network: downsampling inputs] ###
    # Entry block
    x = layers.Conv2D(32, 3, strides=2, padding="same")(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Activation("relu")(x)

    previous_block_activation = x  # Set aside residual

    # Blocks 1, 2, 3 are identical apart from the feature depth.
    for filters in [64, 128, 256]:
        x = layers.Activation("relu")(x)
        x = layers.SeparableConv2D(filters, 3, padding="same")(x)
        x = layers.BatchNormalization()(x)

        x = layers.Activation("relu")(x)
        x = layers.SeparableConv2D(filters, 3, padding="same")(x)
        x = layers.BatchNormalization()(x)

        x = layers.MaxPooling2D(3, strides=2, padding="same")(x)

        # Project residual
        residual = layers.Conv2D(filters, 1, strides=2, padding="same")(
            previous_block_activation
        )
        x = layers.add([x, residual])  # Add back residual
        previous_block_activation = x  # Set aside next residual

    ### [Second half of the network: upsampling inputs] ###

    for filters in [256, 128, 64, 32]:
        x = layers.Activation("relu")(x)
        x = layers.Conv2DTranspose(filters, 3, padding="same")(x)
        x = layers.BatchNormalization()(x)

        x = layers.Activation("relu")(x)
        x = layers.Conv2DTranspose(filters, 3, padding="same")(x)
        x = layers.BatchNormalization()(x)

        x = layers.UpSampling2D(2)(x)

        # Project residual
        residual = layers.UpSampling2D(2)(previous_block_activation)
        residual = layers.Conv2D(filters, 1, padding="same")(residual)
        x = layers.add([x, residual])  # Add back residual
        previous_block_activation = x  # Set aside next residual

    # Add a per-pixel classification layer
    outputs = layers.Conv2D(num_classes, 3, activation="softmax", padding="same")(x)

    # Define the model
    model = keras.Model(inputs, outputs)
    return model


# Free up RAM in case the model definition cells were run multiple times
keras.backend.clear_session()

# Instantiate data Sequences for each split
IMG_SIZE = (128, 128)
train_gen = TrainingData(1, IMG_SIZE, train_input_img_paths, train_target_img_paths)
val_gen   = TrainingData(1, IMG_SIZE, val_input_img_paths, val_target_img_paths)

num_classes = 2
batch_size = 32

# model for current phase only
model = get_model(IMG_SIZE, 2)
model.compile(optimizer="rmsprop", loss="sparse_categorical_crossentropy")

callbacks = [
    keras.callbacks.ModelCheckpoint("bdello_segmentation.h5", save_best_only=True, monitor='val_loss', mode='min'),
    keras.callbacks.EarlyStopping(monitor='val_loss', patience=20, verbose=0, mode='min'),
]

# model fortimelapse 
# give a 3D image of [last_frame_segmentation, current_frame, current_phase_segmentation]
# current_phase_segmentation -> output from previous model
model_tl = get_model(IMG_SIZE, 2)
model_tl.compile(optimizer="rmsprop", loss="sparse_categorical_crossentropy")

# model fortimelapse alt
# give a 2D image of [last_frame_segmentation, current_frame]
model_tl2 = get_model((128,128,2), 2)
model_tl2.compile(optimizer="rmsprop", loss="sparse_categorical_crossentropy")



# Train the model, doing validation at the end of each epoch.
epochs = 100
model.fit(train_gen, epochs=epochs, validation_data=val_gen, callbacks=callbacks)



val_gen = TrainingData(1, IMG_SIZE, val_input_img_paths, val_target_img_paths)
val_preds = model.predict(val_gen)

model.save(WORK_PATH)


def display_mask(i):
    """Quick utility to display a model's prediction."""
    mask = np.argmax(val_preds[i], axis=-1)
    mask = np.expand_dims(mask, axis=-1)
    img = PIL.ImageOps.autocontrast(keras.preprocessing.image.array_to_img(mask))
    return plt.imshow(img, plt.get_cmap("binary_r"))


i = 1
def showresult(i):
    plt.subplot(1,3,1)

    img = imread(val_input_img_paths[i])
    phase = img[1,:,:]
    binary = img[0,:,:]
    plt.imshow(phase, cmap= plt.get_cmap("magma_r"))
    plt.axis('off')
    plt.title('Phase')
    plt.subplot(1,3,2)
    # img = Image.open(val_fluo_img_paths[i])
    # plt.imshow(np.asarray(img), plt.get_cmap("magma"))
    # plt.axis('off')
    # plt.title('Signal')
    # plt.subplot(1,4,3)
    img = Image.open(val_target_img_paths[i])
    plt.imshow(binary, plt.get_cmap("binary_r"))
    plt.axis('off')
    plt.title('Expected\nsegmentation')
    plt.subplot(1,3,3)
    display_mask(i)
    plt.axis('off')
    plt.title('CNN\nsegmentation')
    plt.show()

for i in range(50,100):
    showresult(i)

showresult(8)

# Display results for validation image #10


# Display input image

val_target_img_paths[66]
# Display ground-truth target mask

plt.imshow(img)

# Display mask predicted by our model
display_mask(i)  # Note that the model only sees inputs at 150x150.



model = keras.models.load_model(Path(Path.home()/"CNN"/"Bdello_UNet"))
model.summary()
