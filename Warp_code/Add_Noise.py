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
synthetic_data_fldr = seq_tag + "_" + frame_tag + '_SyntheticTextures/'
synthetic_data = 'Homog_MaxAct0.1_VFA_G16x16/'
synthetic_fr_fldr = synthetic_data_fldr + synthetic_data
synthetic_all_frames = [file for file in os.listdir(synthetic_fr_fldr) if file.endswith('.tif')]
synthetic_all_frames = sorted(synthetic_all_frames, key=lambda f: int(re.sub('\D', '', f)))

# Read first frame
first_frame = io.imread(synthetic_fr_fldr + synthetic_all_frames[0])
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

all_raw_frames = [io.imread(synthetic_fr_fldr + synthetic_all_frames[ii]) for ii in range(len(synthetic_all_frames))]
all_frames_uint8 = fcns.proces_uint16_uint8(all_raw_frames, mask)

# Define list of magnitude ratios and octaves for Perlin noise
ratio_lst = np.linspace(1,10,10)
ratio_lst = ratio_lst*0.01
octv_lst = np.linspace(4,80,20)

# Add Perlin noise 
for mag_r in ratio_lst:
    for octv in octv_lst:
        # Create empty list to save noisy frames for each magnitude ratio and octave setting 
        all_noisy_frames = []
        
        # Create folder to save noisy frames to folder
        noisy_frame_folder = os.path.join('Noisy_' + seq_tag + "_" + frame_tag + '_SyntheticTextures/','Homog_MaxAct0.1_VFA_G16x16_MagR{0}_Oct{1}/'.format(mag_r,octv),"")
        if not os.path.exists(noisy_frame_folder):
            os.makedirs(noisy_frame_folder)

        img_files = glob.glob(noisy_frame_folder + '/*.tif')
        for ff in img_files:
            os.remove(ff)
        
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

            imageio.imwrite(noisy_frame_folder + 'Noisy_Frame%04d.tif'%(fp),noisy_frame)
        
        # Save noisy synthetic texture as gif
        video_name = 'Noisy_Synthetic_Perlin_MagR{0}_Oct{1}.gif'.format(mag_r,octv)    
        writer = imageio.get_writer(noisy_frame_folder + video_name, fps=10)
        for noisy_img in all_noisy_frames:
            writer.append_data(noisy_img)
        writer.close()
        