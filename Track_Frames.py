# =============================================================================
# Feature tracking in region of interest 
# =============================================================================
import matplotlib.pyplot as plt
from skimage import io
import numpy as np
import fcns
import os
import re

# =============================================================================
# Track real and synthetic data: 
#       Compare drift in synthetic data with real estimated drift
# =============================================================================

# Path to real frame images
real_frames_path = 'Frames/' 
seq_tag = ['D1T3_Before','D1T4_Before','D2T4_Before','D2T5_After','D2T6_Before']

# Load FEA results 
data_folder = 'FEA_Results_Homog_MaxAct0.1_VFA/'
disp_all_files = [file for file in os.listdir(data_folder) if file.startswith('disp')]
disp_all_files = sorted(disp_all_files, key=lambda f: int(re.sub('\D', '', f)))

""" Track real data: Loop over all available data """
###################################
track_folder = 'Tracking_Results/' 
if not os.path.exists(track_folder):
    os.mkdir(track_folder)
###################################
all_real_slope = []
for tag in seq_tag:  
    # Load synthetic mask image
    mask_folder = 'Masks/' 
    mask = io.imread(mask_folder + 'Synthetic_Mask_' + tag +'.tif')
    
    # Load tissue mask image
    mask_folder = 'Masks/' 
    tissue_mask = io.imread(mask_folder + tag +'_Mask.tif')
    image_height, image_width = tissue_mask.shape
    tiss_ind = np.argwhere(tissue_mask > 0)
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
    pts_indx = np.argwhere((xy.T[0] > 0.8) & (xy.T[0] < 1.2)  
                           & (xy.T[1] > 0.17) & (xy.T[1] < 0.57))
    pts_indx = np.ravel(pts_indx)

    # Find mean absolute displacement 
    all_gt_disp = [np.loadtxt(data_folder + disp_all_files[nn], dtype=np.float32) for nn in range(len(disp_all_files))]
    all_gt_abs_disp = [np.hypot(all_gt_disp[ff][0,pts_indx]/(scale_factor*conv_factor),all_gt_disp[ff][1,pts_indx]/(scale_factor*conv_factor)) for ff in range(len(all_gt_disp))]
    mean_disp_mid_gt = np.mean(all_gt_abs_disp,axis=1)
    
    real_frames_folder = real_frames_path + 'Frames_' + tag + '/'
    real_frame_imgs = [img for img in os.listdir(real_frames_folder) if img.startswith('Frame')]
    real_frame_imgs = sorted(real_frame_imgs, key=lambda f: int(re.sub('\D', '', f)))
    all_real_frames = [io.imread(real_frames_folder + real_frame_imgs[ii]) for ii in range(len(real_frame_imgs))]
    all_real_frames_uint8 = fcns.proces_uint16_uint8(all_real_frames, mask)
    
    real_mean_abs_disp_track, real_peaks,_, real_mid_valleys = fcns.feature_track_roi(all_real_frames_uint8, mask)
    
    real_mean_abs_disp_valleys = [real_mean_abs_disp_track[vv] for vv in real_mid_valleys]
    real_slope,_ = fcns.find_slope(real_mid_valleys,real_mean_abs_disp_valleys)
    all_real_slope.append(real_slope)

    ###############################################
    # Save results for real  data
    file_name_real = 'Results_Real_Data_' + tag + '.npz'
    np.savez(track_folder + file_name_real, magnitude_ratio=0, octave=0, mean_absolute_displacement=real_mean_abs_disp_track, peaks=real_peaks, valleys=real_mid_valleys, slope=real_slope)
    ##############################################

mean_real_slope = np.mean(all_real_slope)
std_dev_real_slope = np.std(all_real_slope)
##############################################
# Save average results for real  data
file_name_real = 'Average_Results_Real_Data.txt'
header = 'Mean_Slope\nStandard_Deviation'
results =[mean_real_slope,std_dev_real_slope]
np.savetxt(track_folder + file_name_real, results, delimiter=',', fmt='%.6f', header=header)
##############################################

""" Track synthetic data with no noise """
seq_tag_synth = 'D1T3_Before'
frame_tag = 'Frame0000'

syn_fldr = seq_tag_synth + "_" + frame_tag + '_SyntheticTextures/Homog_MaxAct0.1_VFA_G16x16_Square/'
syn_frames_files = [file for file in os.listdir(syn_fldr) if file.endswith('.tif')]
syn_frames_files = sorted(syn_frames_files, key=lambda f: (re.sub('\D', '', f)))

all_syn_frames =  [io.imread(syn_fldr + syn_frames_files[ii]) for ii in range(len(syn_frames_files))]
all_syn_frames_uint8 = fcns.proces_uint16_uint8(all_syn_frames,mask)
syn_mean_abs_disp_track, syn_peaks,_, syn_mid_valleys = fcns.feature_track_roi(all_syn_frames_uint8, mask)

syn_mean_abs_disp_valleys = [syn_mean_abs_disp_track[vv] for vv in syn_mid_valleys]
syn_slope,_ = fcns.find_slope(syn_mid_valleys,syn_mean_abs_disp_valleys)

###############################################
# Save results for synthetic data with no noise
file_name_syn = 'Results_Synthetic_Data_NoNoise_' + seq_tag_synth + '_' + frame_tag + '.npz'
np.savez(track_folder + file_name_syn, magnitude_ratio=0, octave=0, mean_absolute_displacement=syn_mean_abs_disp_track, peaks=syn_peaks, valleys=syn_mid_valleys, slope=syn_slope)
###############################################

""" Track noisy synthetic data """
# Load noisy synthetic frame images 
path_noisy_data = 'Noisy_' + seq_tag_synth + "_" + frame_tag + '_SyntheticTextures/'
all_noisy_syn_frames_folders = [folder for folder in os.listdir(path_noisy_data) if folder.startswith('Homog')]
all_noisy_syn_frames_folders = sorted(all_noisy_syn_frames_folders, key=lambda f: (re.sub('\D', '', f)))

# Loop over the all noisy synthetic data folders and track the set of frames
# Define empty lists to store values of interest from all noisy synthetic data folders
all_nsyn_mean_abs_disp_track = []
all_nsyn_mid_valleys = []
all_nsyn_peaks = []
all_perc_error_slope = []
   
for nsyn_fldr in all_noisy_syn_frames_folders:
    nsyn_frame_imgs = [img for img in os.listdir(path_noisy_data + nsyn_fldr) if img.endswith('tif')]
    nsyn_frame_imgs = sorted(nsyn_frame_imgs, key=lambda f: int(re.sub('\D', '', f)))
    
    all_nsyn_frames =  [io.imread(path_noisy_data + nsyn_fldr +'/'+ nsyn_frame_imgs[ii]) for ii in range(len(nsyn_frame_imgs))]

    nsyn_mean_abs_disp_track, nsyn_peaks,_, nsyn_mid_valleys = fcns.feature_track_roi(all_nsyn_frames, mask)
   
    all_nsyn_mean_abs_disp_track.append(nsyn_mean_abs_disp_track)
    all_nsyn_mid_valleys.append(nsyn_mid_valleys)
    all_nsyn_peaks.append(nsyn_peaks)
    
    nsyn_mean_abs_disp_valleys = [nsyn_mean_abs_disp_track[vv] for vv in nsyn_mid_valleys]
    nsyn_slope,_ = fcns.find_slope(nsyn_mid_valleys,nsyn_mean_abs_disp_valleys)
    
    nsyn_slope_error = (abs(nsyn_slope - mean_real_slope)/mean_real_slope)*100
    all_perc_error_slope.append(nsyn_slope_error)
    
    
    ###############################################
    # Save results for noisy synthetic data
    MR, OCT = fcns.get_options_fl(nsyn_fldr)
    file_name_nsyn = 'Results_Noisy_Synthetic_Data_' + seq_tag_synth + '_' + frame_tag+ '_MR{0}_OCTV{1}.npz'.format(float(MR),int(OCT))
    np.savez(track_folder + file_name_nsyn, magnitude_ratio=float(MR), octave=int(OCT), mean_absolute_displacement=nsyn_mean_abs_disp_track, peaks=nsyn_peaks, valleys=nsyn_mid_valleys, slope=nsyn_slope, slope_error=nsyn_slope_error)
    ###############################################
    
# Find the noisy synthetic data with minimum drift error
min_error_indx = np.argmin(all_perc_error_slope)
min_error_MR, min_error_OCT = fcns.get_options_fl(all_noisy_syn_frames_folders[min_error_indx])

min_nsyn_mean_abs_disp_track_array = np.array(all_nsyn_mean_abs_disp_track[min_error_indx])
min_nsyn_peaks_array = np.array(all_nsyn_peaks[min_error_indx])
min_nsyn_valleys_array = np.array(all_nsyn_mid_valleys[min_error_indx])

# Plot mean absolute displacement corresponding to noisy data with minimum error
plt.figure(figsize=(6, 3.5))
plt.title('Mean Absolute Displacement - Lowest Drift Error')
plt.plot(mean_disp_mid_gt, 'k:',label='Ground Truth FEA')
plt.plot(min_nsyn_mean_abs_disp_track_array, 'teal' ,label='Noisy Synthetic Data')
plt.scatter(min_nsyn_peaks_array,min_nsyn_mean_abs_disp_track_array[min_nsyn_peaks_array],c='firebrick',s=40,marker='o',label='Peaks')
plt.scatter(min_nsyn_valleys_array,min_nsyn_mean_abs_disp_track_array[min_nsyn_valleys_array],c='salmon',s=40,marker='o',label='Valleys')
plt.xlabel('Frame')
plt.ylabel('Mean Displacement')
plt.legend(loc=[0.05,-0.3],ncol=4, prop={'size':7})
plt.show()