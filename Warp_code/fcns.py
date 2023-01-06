"""This document contains all the user-defined functions used for 
generating synthetic data and adding Perlin noise"""

from perlin_noise import PerlinNoise
from skimage import io
import numpy as np
import imageio
import glob
import re
import os

# Extract individual frames from .tif movie files 
def tiff_to_imgs(movie_path,movie_tag):
    """
    Input: movie_path -- Path to the location of the .tif movie files 
           movie_tag -- Reference to the movie filename (for example 'D1T3_Before')
    Output: Tif images of individual movie frames saved in 'save_path'
            save_path -- Created path to save individual movie frames
    """
    save_path = os.path.join('Frames/', 'Frames_' + movie_tag, "")
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    files = glob.glob(save_path + '/*.tif')
    for f in files:
        os.remove(f)
        
    tiff_imgs = [img for img in os.listdir(movie_path) if img.startswith(movie_tag) and img.endswith('.tif')]
    tiff_imgs = sorted(tiff_imgs, key=lambda f: int(re.sub('\D', '', f)))

    count = 0
    for image in tiff_imgs:
        # read image stack
        im = io.imread(movie_path + image)
        nmbr_frames = im.shape[0]
    
        for frame_kk in range(0,nmbr_frames):
            # Saving the images as 16-bit
            imageio.imwrite(save_path + '/Frame%04d.tif'%(frame_kk + count),im[frame_kk])
        count += nmbr_frames
    return save_path

# Calculate deformed 2D coordinates from projective transformation
def find_def_coord_proj_tform(ref_pts,transformation):
    """
    Input: ref_pts -- 2xn matrix of points [x,y] in reference configuration
           transformation -- Estimated projective transformation from skimage
    Output: 2xn matrix of points [x,y] in transformed (deformed) configuration 
    """
    # Write reference coordinates in homogeneous representation
    zz = np.ones((1,np.shape(ref_pts)[1]))
    ref_pts_hom = np.concatenate([ref_pts,zz])
    # Get transformation parameters
    H = transformation.params
    # Find deformed coordinates
    pts_def_scaled = H.dot(ref_pts_hom)
    elation = pts_def_scaled[2,:]
    pts_def = pts_def_scaled/elation 
    return pts_def[0:2,:]

# Process uint16 images --> uint8 images: Bound maximum and minimum pixel intensity according to mask region
def proces_uint16_uint8(frame_list, mask_region): 
    """
    Input: frame_list -- A list of uint16 mxn matrices representing frame images to be processed
           mask_region -- A mxn binary matrix with zeros marking the background and 255 the 
                          region of interest 
                          **If the whole image is a region of interest then the mask would be 
                          a mxn matrix of 255
    Output: A list of uint8 mxn matrices representing processed images
    """
    max_frame = []
    min_frame = []
    
    all_bounded_frames = []
    all_frames_uint8 = []
    
    for frame_indx in range(len(frame_list)):
        frame = frame_list[frame_indx]
        
        masked_region_indx = np.argwhere(mask_region > 0)
        start_crop_y_x = np.min(masked_region_indx,axis=0)
        end_crop_y_x = np.max(masked_region_indx,axis=0)
        cropped_roi  = frame[start_crop_y_x[0]:end_crop_y_x[0]+1,start_crop_y_x[1]:end_crop_y_x[1]+1]  
        cropped_roi_1d = cropped_roi.ravel()
        # Get non-zero pixel values 
        non_zero_indx = np.argwhere(cropped_roi_1d > 0)
        non_zero_frame_values = cropped_roi_1d[non_zero_indx]
        max_frame.append(max(non_zero_frame_values))
        min_frame.append(min(non_zero_frame_values))
    
    max_val = max(max_frame)
    min_val = min(min_frame)
    
    for frame_indx in range(len(frame_list)):
        frame = frame_list[frame_indx]
        bounded_frame = frame.copy()
        bounded_frame = bounded_frame.clip(min = min_val[0] , max = max_val[0])
        reduced_frame = bounded_frame/max_val[0]*255
        reduced_frame = reduced_frame.astype('uint8')
        all_bounded_frames.append(bounded_frame)
        all_frames_uint8.append(reduced_frame)
    
    return all_frames_uint8

# Add Perlin noise to input image
def generate_perlin_noise(img,octaves,seed,mag_ratio):
    """
    Input: img -- A mxn matrix representing the frame image to which noise is to be added
           octaves -- A positive float to assign the number of sub rectangles in each [0, 1] range
           seed -- A positive integer to seed the random integer generator  
           mag_ratio -- A float in the range [0, 1] to assign the intesity of the added noise with 
                        respect to the original image
    Output: A mxn matrix representing the original frame image with added noise 
    """
    noise = PerlinNoise(octaves,seed)
    pix1, pix0 = img.shape
    pnoise = [[noise([i/pix0, j/pix1]) for j in range(pix0)] for i in range(pix1)]
    # make perlin noise from range 0-1
    pnoise = (pnoise - np.min(pnoise)) / (np.max(pnoise) - np.min(pnoise))
    max_image = np.max(img)
    img_rgn_with_noise = img + pnoise * max_image * mag_ratio 
    return img_rgn_with_noise
          