# =============================================================================
# Generate synthetic data based on FEA results
# =============================================================================
from skimage.transform import warp
import matplotlib.pyplot as plt
from skimage import transform
from functools import reduce
from skimage import io
from PIL import Image
import numpy as np
import imageio
import time
import glob
import fcns
import os
import re

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

# Load FEA results 
data_folder = 'FEA_Results_Homog_MaxAct0.1_VFA/'
disp_all_files = [file for file in os.listdir(data_folder) if file.startswith('disp')]
disp_all_files = sorted(disp_all_files, key=lambda f: int(re.sub('\D', '', f)))

# Get the 2D coordinates of points in the middle region (ROI: Square mask) in reference frame 0 (FEA model)
pts_indx = np.argwhere((xy.T[0] > 0.8) & (xy.T[0] < 1.2)  
                       & (xy.T[1] > 0.17) & (xy.T[1] < 0.57))
pts_indx = np.ravel(pts_indx)
pts_0_mid = xy.T[:,pts_indx]

pts_0_mid = pts_0_mid/(scale_factor*conv_factor)

for ss in range(len(shift_factor_X_Y)):
    pts_0_mid.T[:,ss] = pts_0_mid.T[:,ss] + shift_factor_X_Y[ss]
    
# Crop texture image to middle region only (keep dimensions of original image)
strt_crop_x = int(min(pts_0_mid[0]))
end_crop_x = int(max(pts_0_mid[0]))
strt_crop_y = int(min(pts_0_mid[1]))
end_crop_y = int(max(pts_0_mid[1]))

mask = np.zeros([image_height, image_width])
mask[strt_crop_y:end_crop_y, strt_crop_x:end_crop_x] = 255
mask = mask.astype('uint8')
imageio.imwrite(mask_folder + 'Synthetic_Mask_' + seq_tag +'.tif',mask)

# To add border padding when needed
pad = 6 
strt_crop_x_pad = strt_crop_x - pad
end_crop_x_pad = end_crop_x + pad
strt_crop_y_pad = strt_crop_y - pad
end_crop_y_pad = end_crop_y + pad

middle_region = np.zeros([image_height, image_width], dtype=np.uint16)
middle_region[strt_crop_y:end_crop_y,strt_crop_x:end_crop_x] = texture_image[strt_crop_y:end_crop_y,strt_crop_x:end_crop_x]

# List of grid sizes
grid_size = [[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9],[10,10],[11,11],[12,12],[13,13],[14,14],[15,15],[16,16]]

# Overlap tolerance
tol = 2

# Number of frames in FEA simualtion
nmbr_steps = len(disp_all_files)

# Create empty lists to collect results from all grid sizes
all_grid_size_lst_str = []
all_mean_Error_F_allFrames = []
all_peak_mean_Error_F_allFrames = []
all_peak_Error_F_allFrames = []

# Go over all grid sizes in list
for dd in range(len(grid_size)):
    
    start_time_generate_syn_data = time.perf_counter()
    
    # Divide the Middle Region into subdomains 
    nmbr_div_x = grid_size[dd][0]
    nmbr_div_y = grid_size[dd][1]
    
    end_x = max(pts_0_mid[0,:])
    strt_x = min(pts_0_mid[0,:])
    end_y = max(pts_0_mid[1,:])
    strt_y = min(pts_0_mid[1,:])
    
    all_corner_pts_x = np.linspace(strt_x,end_x,nmbr_div_x+1, dtype=np.float32)
    all_corner_pts_y = np.linspace(strt_y,end_y,nmbr_div_y+1, dtype=np.float32)

    # Create empty list to collect results from each grid size
    grid_mean_Error_F = []
    grid_peak_Error_F = []
    grid_syn_frames = []
    
    grid_size_str = ('{0}x{1}').format(nmbr_div_y,nmbr_div_x)
    all_grid_size_lst_str.append(grid_size_str)

    # Create folder to save warped frames with black border
    warped_frame_folder = os.path.join(seq_tag + "_" + frame_tag + '_SyntheticTextures/', 'Homog_MaxAct0.1_VFA_G{0}/'.format(grid_size_str), "") 
    
    if not os.path.exists(warped_frame_folder):
        os.makedirs(warped_frame_folder)
    
    warped_ff = glob.glob(warped_frame_folder + '/*.tif')
    for wf in warped_ff:
        os.remove(wf)
    
    # Deform texture image by data from each step
    for nn in range (nmbr_steps):
        
        # Load FEA data file (displacement) for frame
        u = np.loadtxt(data_folder + disp_all_files[nn], dtype=np.float32)
        # Scale displacement data
        u = u/(scale_factor*conv_factor) 
        # Create empty lists to collect results from step
        warped_sub_frame = []
        grid_Error_F = []
        
        # Loop over subdomains for each grid size
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
              
                # Get 2D displacement in subdomain
                u_sub = u[0:2,pts_indx_sub]
                
                # Get 2D coordinates of deformed points in subdomain
                pts_def_sub = pts_0_sub + u_sub
            
                # Crop texture image to subdomain
                strt_crop_sub_x = int(x_1)
                end_crop_sub_x  = int(x_2)
                strt_crop_sub_y = int(y_1)
                end_crop_sub_y  = int(y_2)
                
                sub = np.zeros([image_height, image_width], dtype=np.float32)
                sub[strt_crop_sub_y:end_crop_sub_y,strt_crop_sub_x:end_crop_sub_x] = middle_region[strt_crop_sub_y:end_crop_sub_y,strt_crop_sub_x:end_crop_sub_x]
                sub = sub.astype('uint16')
                
                # Find the projective transformation 
                tform = transform.estimate_transform('projective', pts_0_sub.T, pts_def_sub.T)
            
                # Warp each subdomain
                warped_sub = warp(sub, tform.inverse, order = 1, preserve_range=True)
                warped_sub = warped_sub.astype('uint16')
                
                warped_sub_frame.append(warped_sub)
                    
                # Find warped coordinates (deformed points due to transformation)
                pts_def_syn = fcns.find_def_coord_proj_tform(pts_0_sub,tform)
                   
                # Calculate deformation gradient for real and synthetic points
                lamda0 = fcns.find_vectors(pts_0_sub)
                lamda = fcns.find_vectors(pts_def_sub)
                lamda_w = fcns.find_vectors(pts_def_syn)

                F_avg_real_sub = fcns.find_F(lamda0,lamda)
                F_avg_syn_sub = fcns.find_F(lamda0,lamda_w)
                
                # Calculate Error in deformation gradient 
                Error_F_sub = abs(F_avg_real_sub - F_avg_syn_sub)/abs(F_avg_real_sub)*100
                grid_Error_F.append(Error_F_sub.ravel())
                    
            mean_Error_F = np.mean(grid_Error_F, axis=0)
            peak_Error_F = np.max(grid_Error_F, axis=0)
                     
        grid_mean_Error_F.append(mean_Error_F) 
        grid_peak_Error_F.append(peak_Error_F)
        
        #Combine the warped subdomains into a single warped region
        whole_warped_frame = np.zeros([image_height, image_width], dtype=np.uint16)
        for jj in range(len(warped_sub_frame)):
            whole_warped_frame = reduce(np.maximum,[whole_warped_frame,warped_sub_frame[jj]])
            whole_warped_frame = whole_warped_frame.astype('uint16')
        
        grid_syn_frames.append(whole_warped_frame)
    
        # Save warped frames as uint16 images (all frames)
        imageio.imwrite(warped_frame_folder + 'Warped_Frame%04d.tif'%(nn),whole_warped_frame)

    run_time = time.perf_counter() - start_time_generate_syn_data  
    
    mean_Error_F_allSub_allFrames = np.mean(grid_mean_Error_F,axis=0)
    all_mean_Error_F_allFrames.append(mean_Error_F_allSub_allFrames)
    
    peak_mean_Error_F_allSub_allFrames = np.max(grid_mean_Error_F,axis=0)
    all_peak_mean_Error_F_allFrames.append(peak_mean_Error_F_allSub_allFrames)
    
    peak_Error_F_allSub_allFrames = np.max(grid_peak_Error_F,axis=0)
    all_peak_Error_F_allFrames.append(peak_Error_F_allSub_allFrames)
    
    # Save synthetic data as gif (all frames)
    video_name = warped_frame_folder + 'Synthetic_Frames_G{0}.gif'.format(grid_size_str)
    writer = imageio.get_writer(video_name, fps=10)
    for image in grid_syn_frames:
        writer.append_data(image)
    writer.close()

# =============================================================================
# # Plot Errors
# =============================================================================
# Plot Peak-Mean Error
plt.figure(figsize=(10, 5), dpi=300)
x = range(len(grid_size))
plt.xticks(x, all_grid_size_lst_str)
y_p = np.array(all_peak_mean_Error_F_allFrames)
plt.plot(x, y_p[:,0],marker ='s',label='$F_{xx}$')
plt.plot(x, y_p[:,1],marker ='d',label=r'$F_{xy}$')
plt.plot(x, y_p[:,2],marker ='.',label=r'$F_{yx}$')
plt.plot(x, y_p[:,3],marker ='x',label=r'$F_{yy}$')
plt.grid('on')
plt.xlabel('Grid Size')
plt.ylabel('Percentage Peak-Mean Error')
plt.legend()
plt.show()


# Plot Peak Error
plt.figure(figsize=(10, 5), dpi=300)
x = range(len(grid_size))
plt.xticks(x, all_grid_size_lst_str)
y_p = np.array(all_peak_Error_F_allFrames)
plt.plot(x, y_p[:,0],marker ='s',label='$F_{xx}$')
plt.plot(x, y_p[:,1],marker ='d',label=r'$F_{xy}$')
plt.plot(x, y_p[:,2],marker ='.',label=r'$F_{yx}$')
plt.plot(x, y_p[:,3],marker ='x',label=r'$F_{yy}$')
plt.grid('on')
plt.xlabel('Grid Size')
plt.ylabel('Percentage Peak Error')
plt.legend()
plt.show()        

