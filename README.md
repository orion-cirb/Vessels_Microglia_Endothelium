# Vessels_Cellpose_Microglia

* **Developed for:** Nicolas
* **Team:** Garel
* **Date:** July 2024
* **Software:** Fiji

### Images description

3D images taken with a x20 objective on a confocal microscope.

3 channels:
  1. Vessels (mandatory)
  2. Microglia (optional)
  3. Endothelial nuclei (optional)

With each image can be provided a *.roi* or *.zip* file containing one or multiple ROI(s):
- if no ROI provided, whole image is analyzed
- if ROI *roiName* provided alone, analysis is performed in ROI not translated along slices
- if ROIs *roiName* and *roiName_2* provided, analysis is performed in ROI translated from *roiName* to *roiName_2* position along slices

### Plugin description

* Detect vessels with either:
  * Quantile Based Normalization (optional) + Cellpose + closing + median filtering
  * Quantile Based Normalization (optional) + median filtering + DoG + thresholding + closing + median filtering
* Compute vessels skeleton and provide vessels volume, diameter, length, branches number, junctions numbers, etc.
* If microglia channel provided,
  * Detect cells with median filtering + thresholding + closing
  * Compute distance and volume of contact between each microglial cell and its nearest vessel
* If endothelial nuclei channel provided,
  * Detect nuclei in 2D with median filtering + Omnipose 2D
  * Stitch 2D masks into 3D volume
  * Provide nuclei number


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *vessels* or *vessels2* model 
* **Omnipose** conda environment + *cyto2_omni* pretrained model

### Version history

Version 1 released on July 15, 2024.
