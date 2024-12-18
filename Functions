import pandas as pd
from PIL import Image
import numpy as np
from scipy import stats
import os
from os import listdir
import csv  
import tkinter as tk
from tkinter import filedialog

##Functions

##remove b from a
##a is your image of interest,
##b is the background or autofluorescence you want to remove
##input one frame per image
##nuc = if image is nuclear staining, the output will be a binary image
def rm_image_background(a, b, nuc = False):
    #Convert images to arrays
    im1_image_array = np.array(b)
    im2_image_array = np.array(a)
    #Remove negative values
    im1_array = im1_image_array[im1_array>0]
    im2_array = im2_image_array[im2_array>0]
    #Get the "modes" for the arrays representing the images
    mode_im1_array = stats.mode(im1_array)
    mode_im2_array = stats.mode(im2_array)
    #Ratio between image of interest and background noise
    mode_ratio = mode_im1_array/mode_im2_array
    #Normalizing to background signal
    im2_image_array = im2_image_array*mode_ratio
    #Remove background signal and autofluorescence
    im3_image_array = im2_image_array-im1_image_array
    #Set negative values to 0
    im3_image_array = im3_image_array * (im3_image_array > 0)
    #Remove values from non-overlapping regions
    im3_image_array = im3_image_array * (im1_image_array > 0)
    if nuc:
        test_array = im3_image_array[im3_image_array > 0]
        #Add mirrored data to generate a two-tailed distribution
        full_array = [-x for x in test_array] + test_array
        sd_full_array = np.std(full_array)
        #Extract data above 2 standard deviations from the mean
        im2_binary = im3_image_array>(2*sd_full_array)
        return im2_binary
    else:
        return im3_image_array
    

##bin image, either nuclei segmented image or intensity image
##a = image of interest (one frame/channel)
##bin_width = binsize in pixels
##nuc = if image is nuclei segmented and labeled (for example with bwlabel)
def bin_image(a, bin_width, nuc = False):
    mat = np.array(a)
    bin_size = bin_width
    n_col = np.floor(mat.shape[1]/bin_size)
    n_row = np.floor(mat.shape[0]/bin_size)
    mat_bin = np.zeros(n_row,n_col)
    #Calculates number of cells for each bin
    if nuc:
        tab, cell_id, cell_size = np.unique(mat[mat > 0], return_inverse=True, return_counts=True)
        for i in range(n_row):
            for j in range(n_col):
                k1 = i * bin_size
                k2 = (i + 1) * bin_size
                m1 = j * bin_size
                m2 = (j + 1) * bin_size
                ss = 0
                v = np.unique(mat[k1:k2, m1:m2])
                v1 = v[v!=0]
                if len(v1) > 0:
                    for p in v1:
                        index = np.where(cell_id == p)[0]
                        current_cell_actual_size = cell_size[index]
                        pxl_cell_count_binned = np.sum(mat[k1:k2, m1:m2] == p)
                        ss += pxl_cell_count_binned/current_cell_actual_size[0]
                
                mat_bin[i, j] = ss
    #Calculates mean intensity for each bin        
    else:
        mat = mat * (mat > 0)
        for i in range(n_row):
            for j in range(n_col):
                k1 = i * bin_size
                k2 = (i + 1) * bin_size
                m1 = j * bin_size
                m2 = (j + 1) * bin_size
                mat_bin[i,j] = np.mean(mat[k1:k2, m1:m2])
                
    return mat_bin

##remove background from all filters in two aligned images (pre and post)
##(input three RGB aligned images, DAPI/FITC pre, DAPI/TRITC pre, FITC/TRITC post), 
##the pre and post staining images.
##the blue filter is dapi, and so will be nuclei which will give a binary image.
##Use bwlabel() on the dapi/blue frame/channel to segment/label the output frame containing the binary frame of cells.
def rm_sample_background(a, b, c):
    rgb_img_fitc_pre = np.array(a)
    rgb_img_tritc_pre = np.array(b)
    rgb_img_fitc_trict_post = np.array(c)
    r_img_pre = rgb_img_tritc_pre[:,:,0]
    g_img_pre = rgb_img_fitc_pre[:,:,1]
	b_img_pre = rgb_img_fitc_pre[:,:,2]
    
    r_img_post = rgb_img_fitc_trict_post[:,:,0]
	g_img_post = rgb_img_fitc_trict_post[:,:,1]

    #Normalize and remove background and autofluorescence AND creates binary nuclei image
	dapi = rm_image_background(b_img_pre,g_img_pre,nuc=True)
    #Normalize and remove background and autofluorescence
	fitc = rm_image_background(g_img_post,g_img_pre,nuc=False)	
	tritc = rm_image_background(r_img_post,r_img_pre,nuc=False)
	#Merges the processed images to one RGB image
    final_img = np.stack(tritc, fitc, dapi, axis = -1)	
	return final_img 


##Function for autobinning of segmented nuclearstained images
##Takes a segmented/labeled image and finds optimal binning size
##by sequential binning and finding a maxima
##Output is the optimal binning size
def autobinning(img, name):
    #img = np.array(x)
    tab, cell_id, cell_size = np.unique(mat[mat > 0], return_inverse=True, return_counts=True)
    # Store information about the autobinning in a .csv-file
    cell_metrics = pd.DataFrame({'cellID': unique_cells, 'cellsize': cellsize})
    csv_filename = f"{name}_cellmetrics.csv"
    cell_metrics.to_csv(csv_filename, index=False)
    
    #Calculate celldiameter from mean cell area (pixels)
    cell_diameter = 2*np.sqrt(np.mean(cell_size[:-1])/np.pi)
    increment = np.ceiling(cell_diameter)
    #Variables and vectors to create output data frame and keep track on bin metrices
	bin_size = 0 		
	current_bincount = 0
	previous_bincount = 0
    #Maximum iterations before terminating the autobinning
	max_iteration = 20								
	bin_over_1 = []
	bin_over_0 = []
	binsizes = []
	tot_bins = []
	counter=0
	previous_slope = 0
	current_slope = 1
	#The automatic binning, keeps iterating until the number of bins containing 
    #multiple cells declines or if the incline is too small
    #while (current_slope > 0 && counter<max_iteration){
    while current_slope > 0 and counter < max_iteration:
        counter += 1							
		print("autobinning iteration number:")
		print(counter)
        #Increase bin size one celldiameter per iteration
		bin_size = bin_size+increment 						
		print("bin width")
		print(bin_size)

		n_col=np.floor(img.shape[1]/bin_size)
		n_row=np.floor(img.shape[2]/bin_size)
        #Create empty matrix for storing the output from binning of the image of current iteration
		mat_bin = np.zeros(n_row,n_col)			
		mat_bin = bin_image(img,bin_size,nuc=True)

		previous_bincount = current_bincount
		previous_slope = current_slope

        #Appends bin size to list of bin sizes
        binsizes.append(bin_size)
        #Calculates the number of bins containing more than 1 cell
		current_bincount = np.sum(mat_bin>1) 	
        #Appends the bin count of multicellular bins to vector of bin counts		
		bin_over_1 = bin_over_1.append(current_bincount)		
        #Calculate total number of bins containing cells
		bin_count_tot = np.sum(mat_bin>0)
        #Appends the total number of bins containing cells to vector of total bins containing cells counts 				
		bin_over_0 = bin_over_0.append(bin_count_tot)		
		#Calculates total number of bins
		tot_bin = np.sum(mat_bin) 						
        #Appends total number of bins to vector of total number of bins
		tot_bins = tot_bins.append(tot_bin) 					
        #Calculate current slope
		current_slope = (current_bincount-previous_bincount)/increment 	
		print("current slope")
		print(current_slope)
	}
	#Summarize all info into a data frame
    df_autobin = pd.DataFrame({
    "bin_over_1": bin_over_1,
    "bin_over_0": bin_over_0,
    "tot_bins": tot_bins,
    "binsizes": binsizes
    })    
    # Returns data fram with information about the bin counts
	return df_autobin 								

##Bin an RGB image, bins the channels/frames seperately then returns a binned RGB image.
##input:
##a = RGB image, representing TRITC/FITC/DAPI, in that order
##bin_width = the bin size in pixels
def bin_rgb_image(a, bin_width):
    #store image as imagedata
	rgb_image = np.array(a)			
	#Extract red filter from image
	r_img = rgb_image[,,0]		
    #Extract green filter from image
	g_img = rgb_image[,,1]			
    #Extract blue filter from image
	b_img = rgb_image[,,2]	
    #Applies the binning function on each channel, DAPI channel
	bin_dapi = bin_image(b_img,bin_width,nuc=True)	
	bin_fitc = bin_image(g_img,bin_width,nuc=False)  #FITC channel
	bin_tritc = bin_image(r_img,bin_width,nuc=False) #TRITC channel
	#Creates and RGB image from the three binned frames
    final_bin_img = np.stack(bin_tritc, bin_fitc, bin_dapi, axis = -1)	
    return final_bin_img


def choose_images:
    root = tk.Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    files = filedialog.askopenfilename(multiple=True) 
    ##%gui tk
    var = root.tk.splitlist(files)
    filepath = []
    for f in var:
        filepath.append(f)
    return filepath

def get_images(filepaths):
    images = []
    for filepath in filepaths:
        images.append(Image.open(filepath))
    return images



def order_images(images):
    channels = ['tritc', 'fitc', 'dapi']
    variants = ['pre', 'post']
    combinations = product(variant, channels)
    number_of_images = len(images)
    indeces = []
    
    for image in images:
        index = 0
        for combination in combinations:  
            if all(x in image.filename for x in combination):
                indeces.append(index)
                break
            index++
    images = [images[i] for i in indeces]
    return images

def merge_images(images):
    merged_image_pre = np.stack(image[0], image[1], image[2], axis = -1)
    merged_image_post = np.stack(image[3], image[4], image[5], axis = -1)	
	return [merged_image_pre, merged_image_post]

def correct_sample_images(images):
    dapi_fitc_pre = np.stack(images[0][,,1], images[0][,,1], images[0][,,2], axis = -1)
    dapi_tritc_pre = np.stack(images[0][,,2], images[0][,,2], images[0][,,2], axis = -1)
    fitc_tritc_post = np.stack(images[1][,,0], images[1][,,1], images[1][,,2], axis = -1)

    correcteded_sample_image = rm_sample_background(dapi_fitc_pre, dapi_tritc_pre, fitc_tritc_post)
    return correcteded_sample_image 


def correct_samples(x):

    root = tk.Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    folder_selected = filedialog.askdirectory()
    print(folder_selected)
    dir_paths = []

    for directoryname in os.listdir(folder_selected):
        directorypath = os.path.join(folder_selected, directoryname)
        if os.path.isdir(directorypath):
            dir_paths.append(directorypath)
            for directory in dir_paths:
                file_paths []
                for filename in os.listdir(directory):
                    filepath = os.path.join(directory, filename)
                    if os.path.isfile(filepath):
                        file_paths.append(filepath)
                        print(filepath)
                        corrected_image = correct_sample_images(merge_images(order_images(get_images(x))))
                        corrected_image.save(f"{os.path.dirname(directory)}/{corrected_image.name}.jpeg")


##From a folder containing samples, each sample having their own folder,
##where each sample needs to have all of its channels in that specific folder.
##all images need to be named by channel with suffix pre or post
##(representing before or after staining with antibody)
##ex. samplename_dapi_pre.jpg, samplename_dapi_post.jpg, 
##samplename_fitc_pre.jpg and samplename_fitc_post.jpg
##the folder needs to contain the DAPI channel and at least one of FITC or TRITC
##
##- not yet - (if you use one marker only, 
##the code will anyway just copy the FITC values into the TRITC or vice verca)
## ex...
## //path/folder_containing_individual_folders_for_each_sample
## //path/folder_containing_individual_folders_for_each_sample/folder_timepoint_1
## //path/folder_containing_individual_folders_for_each_sample/folder_timepoint_1/samplename_dapi_pre.jpg, samplename_dapi_post.jpg, samplename_fitc_pre.jpg and samplename_fitc_post.jpg

def normalize_with_autobinning_samples():
    root = tk.Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    folder_selected = filedialog.askdirectory()
    print(folder_selected)
	# Keetp track of the maximum value for each channel across all samples
    max_values = {"tritc": 0, "fitc": 0, "dapi": 0}
	
    for filename in os.listdir(folder_selected):
        filepath = os.path.join(folder_selected, filename)
        if os.path.isfile(filepath):
            image = Image.open(filepath).convert("RGB")  # Ensure image is in RGB format
            img_to_process = image[,,2]
            autobinned_data = autobinning(img_to_process,image.name) 
            

            # Save autobinned data to a csv.file
            filename = f"{os.path.basename(image.filename)}_autobinned.csv"
            csv_filepath = os.path.join(folder_selected, filename)
            csv_filename = f"{current_image_name}_autobin.csv"
            df_autobin.to_csv(filename, index=False) 


            optimal_bin_width = df_autobin.loc[df_autobin[bin_over_1].idmax(), binsizes]
            binned_RGB = binRGBImage(np.array(image), optimal_bin_width)
            
            #Store metric data 
            metrics_data = {
                "NCOL": binned_RGB.shape[1],
                "NROW": binned_RGB.shape[0],
                "binwidth": optimal_bin_width
            }
            metrics_csv_filename = f"{current_image_name}_metrics.csv"
            metrics_csv_filepath = os.path.join(folder_selected, metrics_csv_filename)
            pd.DataFrame([metrics_data]).to_csv(metrics_csv_filepath, index=False)

            
            binned_RGB_int_per_pixel = binned_RGB/optimal_bin_width
            #Store images as vectors inside a data frame, columns being each channel
            df_normalized_RGB = pd.DataFrame({
                "tritc" : binned_RGB_int_per_pixel[,,0].flatten(),				
                "fitc" : binned_RGB_int_per_pixel[,,1].flatten(),
                "dapi" : binned_RGB_int_per_pixel[,,2].flatten()
            })
            #check maximum value across all samples for each iteration, to update total maximum values
            for key in max_values:
                max_values[key] = max(max_values[key], df_normalized_RGB[key].max())
    #Normalize the adjusted sample images
    for filename in os.listdir(folder_selected):
        filepath = os.path.join(folder_selected, filename)
        if os.path.isfile(filepath):
            image = Image.open(filepath).convert("RGB")  # Ensure image is in RGB format
            img_data = np.array(image)
            
            # Normalize each channel
            normalized_data = np.zeros_like(img_data, dtype=np.float32)
            normalized_data[:, :, 0] = img_data[:, :, 0] / max_values["tritc"]
            normalized_data[:, :, 1] = img_data[:, :, 1] / max_values["fitc"]
            normalized_data[:, :, 2] = img_data[:, :, 2] / max_values["dapi"]

           # Adjust values for image representation
           normalized_data = (normalized_data * 255).clip(0, 255).astype(np.uint8)

           # Convert to an image and save as a new file
           normalized_image = Image.fromarray(normalized_data, mode="RGB")
           normalized_image.save(filepath.replace(".jpg", "_normalized.jpg"))
