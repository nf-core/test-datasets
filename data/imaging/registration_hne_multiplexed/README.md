# modules/nf-core/stainwarpy test data

## Description
This dataset consists of **multiplexed and H&E images of human colon tissue** intended for testing image processing modules in tasks such as multiplexed and H&E image registration/transformation. The images are cropped or simplified versions of real data to allow testing of module functionality.

---

## Contents

- `multiplexed_image_colon.ome.tif`  
  Multiplexed image of a human colon section. Contains two channels `0: DAPI`, `1: IgG1` acquired over different staining cycles. Used to test multi-channel image processing, channel extraction and registration.

- `hne_image_colon.ome.tif`  
  H&E-stained image of the same human colon tissue section. Used for co-registration with corresponding multiplexed image.

- `multiplexed_single_channel_img.ome.tif`  
  A single-channel version of the multiplexed image `0: DAPI`. Useful for co-registration with corresponding H&E image.

- `hne_segmentation_mask.ome.tif`  
  Segmentation mask corresponding to `hne_image_colon.ome.tif`. To be transformed with same transformation for H&E image and useful for overlay, or segmentation-based workflows.

- `feature_based_transformation_map.npy`  
  Pre-computed feature-based transformation map to between the multiplexed and H&E image pair (H&E as moving image). Used to apply the same transformation on segmentation mask.

---

## Source
Derived from real colon tissue images. Files have been cropped and/or simplified to create a small test dataset suitable for module testing. No patient-identifiable information is included.  

---

## Notes
- All image files are in **OME-TIFF** format except the transformation map, which is a NumPy `.npy` array.   
- File sizes are small to allow CI testing.
