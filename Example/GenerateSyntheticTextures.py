# =============================================================================
# Generate synthetic data based on FEA results
# =============================================================================
from skimage.transform import warp
from skimage import transform
from functools import reduce
from skimage import io
from PIL import Image
import numpy as np
import imageio
import fcns
import os

# Path to .tif movie file
seq_fldr = 'Movies/'
# Specify video name tag
seq_tag = 'D1T4_Before'

# Convert movie to individual frames
imgs_fldr = fcns.tiff_to_imgs(seq_fldr,seq_tag)
# Load frame to use as texture image
frame_tag = 'Frame0000'
image_seq = io.imread(imgs_fldr + frame_tag + '.tif')

# Load tissue mask image
mask_folder = 'Masks/' 
tissue_mask = io.imread(mask_folder + seq_tag +'_Mask.tif')
image_height, image_width = tissue_mask.shape

tiss_ind = np.argwhere(tissue_mask > 0)

# Create a tissue texture with transparent background
texture_image = np.zeros([image_height,image_width],dtype='uint16')
transparency_mat = np.zeros([image_height,image_width],dtype='uint16')
for tt in range(len(tiss_ind)):
    ind = tiss_ind[tt]
    texture_image[ind[0],ind[1]] = image_seq[ind[0],ind[1]]
    transparency_mat[ind[0],ind[1]] = 1

# Save the tissue texture as 16-bit grayscale image with transparent background
texture_folder = 'Textures/' 
if not os.path.exists(texture_folder):
    os.mkdir(texture_folder)
array_buffer = texture_image.tobytes()
img = Image.new("I", texture_image.T.shape)
img.frombytes(array_buffer, 'raw', "I;16")
img.save(texture_folder + seq_tag + "_" + frame_tag + "_TissueTexture.png",**{'transparency': 0})

#Create folder to save warped frames 
warped_frame_folder = 'SyntheticData/' 

if not os.path.exists(warped_frame_folder):
    os.makedirs(warped_frame_folder)

# =============================================================================
# Match FEA coordinates with tissue mask
# =============================================================================
# Load tissue slice coordinates (from FEniCS)
xy = np.loadtxt('Tissue_Slice_Coordinates.txt', dtype=np.float32)

# Scale coordinates to match image (originally scaled to have a smaller FE model)
scale_factor = 4
    # Convert um to mm
conv_factor = 1e-3
    # Scale xy coordinates
xy_scale = xy/(scale_factor*conv_factor)
# Align the centers of the mask and the tissue slice coordinates:
    #1 Find center of mask
mask_cntr = np.flip(0.5*np.ptp(tiss_ind,axis=0) + np.min(tiss_ind,axis=0))
    #2 Find center of slice coordinates
slice_cntr = 0.5*np.ptp(xy_scale,axis=0) + np.min(xy_scale,axis=0)

shift_factor_X_Y = mask_cntr - slice_cntr

orig_shape = np.shape(xy)
xy_shifted = np.zeros(orig_shape)

for ss in range(len(shift_factor_X_Y)):
    xy_shifted[:,ss] = xy_scale[:,ss] + shift_factor_X_Y[ss]

xy_compat = xy_shifted.T

# =============================================================================
# Get the 2D coordinates of points in the middle region (ROI: Square mask) in reference frame 0
# =============================================================================
m_ex = 45 # mask extension from center
pts_indx = np.argwhere((xy_compat[0] > (mask_cntr[0] - m_ex)) & (xy_compat[0] < (mask_cntr[0] + m_ex)) 
                        & (xy_compat[1] > (mask_cntr[1] - m_ex)) & (xy_compat[1] < (mask_cntr[1] + m_ex)))
pts_indx = np.ravel(pts_indx)
pts_0_mid = xy_compat[:,pts_indx]
    
# Crop texture image to middle region only (keep dimensions of original image)
end_x = max(pts_0_mid[0,:])
strt_x = min(pts_0_mid[0,:])
end_y = max(pts_0_mid[1,:])
strt_y = min(pts_0_mid[1,:])

mask = np.zeros([image_height, image_width])
mask[int(strt_y):int(end_y),int(strt_x):int(end_x)] = 255
mask = mask.astype('uint8')
imageio.imwrite(mask_folder + 'Synthetic_Mask_' + seq_tag +'.tif',mask)

middle_region = np.zeros([image_height, image_width], dtype=np.uint16)
middle_region[int(strt_y):int(end_y),int(strt_x):int(end_x)] = texture_image[int(strt_y):int(end_y),int(strt_x):int(end_x)]

# =============================================================================
# Load FEA results 
# =============================================================================
data_folder = 'FEA_Results_Homog_MaxAct0.08_VFA_Z_Lite/'
pos_all_files = [file for file in os.listdir(data_folder) if file.startswith('pos')]
pos_all_files = sorted(pos_all_files)

pos_x = np.loadtxt(data_folder + pos_all_files[0])
pos_y = np.loadtxt(data_folder + pos_all_files[1])

nmbr_points,nmbr_steps = pos_x.shape

# Scale FEA points
pos_x_scale = pos_x/(scale_factor*conv_factor) + shift_factor_X_Y[0]
pos_y_scale = pos_y/(scale_factor*conv_factor) + shift_factor_X_Y[1]

# =============================================================================
# Define grid size and warp texture in each grid
# =============================================================================
grid_size = [16,16]

# Overlap tolerance
tol = 2

# Divide the Middle Region into subdomains 
nmbr_div_x = grid_size[0]
nmbr_div_y = grid_size[1]

all_corner_pts_x = np.linspace(strt_x,end_x,nmbr_div_x+1, dtype=np.float32)
all_corner_pts_y = np.linspace(strt_y,end_y,nmbr_div_y+1, dtype=np.float32)

# Create empty list to store synthetic frames for all steps
syn_frames = []
mean_error_pos_step = []
peak_error_pos_step = []

# Deform texture image by data from each step
for nn in range (nmbr_steps):
    
    # position at each step
    pts_def_x = pos_x_scale[:,nn]
    pts_def_y = pos_y_scale[:,nn]
    pts_def = np.array([pts_def_x,pts_def_y])
    
    # Create empty lists to collect results from each step
    warped_sub_frame = []
    mean_error_pos_sub = []
    peak_error_pos_sub = []
    
    # Loop over subdomains 
    for kk in range(len(all_corner_pts_x)-1):
        x_1 = all_corner_pts_x[kk] - tol
        x_2 = all_corner_pts_x[kk+1] + tol
        if kk==0:
            x_1 = all_corner_pts_x[kk]
        if kk==len(all_corner_pts_x)-2:
            x_2 = all_corner_pts_x[kk+1]
        
        for yn in range(len(all_corner_pts_y)-1):
            y_1 = all_corner_pts_y[yn] - tol
            y_2 = all_corner_pts_y[yn+1] + tol
            if yn==0:
                y_1 = all_corner_pts_y[yn]
            if yn==len(all_corner_pts_y)-2:
                y_2 = all_corner_pts_y[yn+1]            
            
            # Get the coordinates of points in subdomain in reference frame 0
            pts_indx_sub = np.argwhere((xy_compat[0] > x_1) & (xy_compat[0] < x_2) & (xy_compat[1] > y_1) & (xy_compat[1] < y_2))
            pts_indx_sub = np.ravel(pts_indx_sub)
            
            pts_0_sub = xy_compat[:,pts_indx_sub]
            
            # Get coordinates of deformed points in subdomain
            pts_def_sub = pts_def[:,pts_indx_sub]
        
            # Crop texture image to subdomain
            sub = np.zeros([image_height, image_width], dtype=np.float32)
            sub[int(y_1):int(y_2),int(x_1):int(x_2)] = middle_region[int(y_1):int(y_2),int(x_1):int(x_2)]
            sub = sub.astype('uint16')
            
            # Find the projective transformation 
            tform = transform.estimate_transform('projective', pts_0_sub.T, pts_def_sub.T)
        
            # Warp each subdomain
            warped_sub = warp(sub, tform.inverse, order = 1, preserve_range=True)
            warped_sub = warped_sub.astype('uint16')
            
            warped_sub_frame.append(warped_sub)
                
            # Find warped coordinates (deformed points due to transformation)
            pts_def_syn = fcns.find_def_coord_proj_tform(pts_0_sub,tform)
            
            pos_sub_error = abs(pts_def_sub - pts_def_syn)/abs(pts_def_sub)*100
            
            mean_sub_pos_error = np.mean(pos_sub_error,axis=1)
            peak_sub_pos_error = np.max(pos_sub_error,axis=1)
        
            
            mean_error_pos_sub.append(mean_sub_pos_error)
            peak_error_pos_sub.append(peak_sub_pos_error)
     
        mean_error_pos = np.mean(mean_error_pos_sub, axis=0)
        peak_error_pos = np.max(peak_error_pos_sub, axis=0)
                
   
    mean_error_pos_step.append(mean_error_pos)
    peak_error_pos_step.append(peak_error_pos)
    
    # Combine the warped subdomains into a single warped region
    whole_warped_frame = np.zeros([image_height, image_width], dtype=np.uint16)
    for jj in range(len(warped_sub_frame)):
        whole_warped_frame = reduce(np.maximum,[whole_warped_frame,warped_sub_frame[jj]])
        whole_warped_frame = whole_warped_frame.astype('uint16')
    
    syn_frames.append(whole_warped_frame)

syn_frames_array = np.array(syn_frames)
imageio.mimwrite(warped_frame_folder + 'SyntheticTexture_14.tif',syn_frames_array,format="tif")
 

all_peak_mean_error_pos = np.max(mean_error_pos,axis=0)
all_peak_error_pos = np.max(peak_error_pos,axis=0)

