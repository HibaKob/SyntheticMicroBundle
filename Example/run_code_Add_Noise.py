# =============================================================================
# Add Perlin noise to generated synthetic data
# =============================================================================
from skimage import io
import numpy as np
import imageio
import glob
import fcns
import os
import re


# Specify video name and frame tags
seq_tag = 'D1T4_Before'
frame_tag = 'Frame0000'

# Load synthetic frames
synthetic_data_fldr = 'SyntheticData/'
synthetic_data = 'SyntheticTexture_14.tif'
synthetic_all_frames = io.imread(synthetic_data_fldr + synthetic_data)


# Read first frame
first_frame = synthetic_all_frames[0]
image_height, image_width = first_frame.shape

# Load synthetic mask image
mask_folder = 'Masks/' 
mask = io.imread(mask_folder + 'Synthetic_Mask_' + seq_tag +'.tif')

# Mask extents
mask_min = np.flip(np.min(np.argwhere(mask>0),axis=0))
mask_max = np.flip(np.max(np.argwhere(mask>0),axis=0))
strt_crop_x = int(mask_min[0])
end_crop_x = int(mask_max[0])+1
strt_crop_y = int(mask_min[1])
end_crop_y = int(mask_max[1])+1

# To add border padding when needed
pad = 6 
strt_crop_x_pad = strt_crop_x - pad
end_crop_x_pad = end_crop_x + pad
strt_crop_y_pad = strt_crop_y - pad
end_crop_y_pad = end_crop_y + pad

all_frames_uint8 = fcns.proces_uint16_uint8(synthetic_all_frames, mask)

# Define magnitude ratio and octave for Perlin noise
#ratio_lst = np.linspace(1,10,10)
mag_r = 0.08
octv = 24

# Add Perlin noise 

# Create empty list to save noisy frames for each magnitude ratio and octave setting 
all_noisy_frames = []

# Create folder to save noisy frames to folder
noisy_frame_folder = 'Noisy_SyntheticData/'
if not os.path.exists(noisy_frame_folder):
    os.makedirs(noisy_frame_folder)

# Initialize seed (for reproducibility)
seed = 200

for fp in range(len(all_frames_uint8)): 
    noisy_frame = np.zeros([image_height, image_width], dtype='uint8') 
    frame_roi = all_frames_uint8[fp][strt_crop_y_pad:end_crop_y_pad,strt_crop_x_pad:end_crop_x_pad]
    noisy_frame_cropped = fcns.generate_perlin_noise(frame_roi,octv,seed,mag_r)
    # Remove noise applied to background (set value below threshold to 0)
    max_noisy = np.max(noisy_frame_cropped)
    threshold = mag_r*max_noisy
    noisy_frame_cropped = np.where(noisy_frame_cropped<threshold,0,noisy_frame_cropped)
    # Shift pixel intensity range between 0 and 255
    noisy_frame_cropped = (noisy_frame_cropped/max_noisy)*255
    noisy_frame_cropped = noisy_frame_cropped.astype('uint8')
    
    noisy_frame[strt_crop_y_pad:end_crop_y_pad,strt_crop_x_pad:end_crop_x_pad] = noisy_frame_cropped
    all_noisy_frames.append(noisy_frame)
    
    # Increment the seed to vary Perlin noise between frames
    seed +=5

all_noisy_frames_array = np.array(all_noisy_frames)
imageio.mimwrite(noisy_frame_folder + synthetic_data[0:-4]+'_MagR{0}_Oct{1}.tif'.format(mag_r,octv) ,all_noisy_frames_array,format="tif")

