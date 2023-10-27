# Vessels_Cellpose_Microglia

* **Developed for:** Nicolas
* **Team:** Garel
* **Date:** October 2023
* **Software:** Fiji

### Images description

3D images taken with a x20 objective on a confocal microscope.

2 channels:
  1. IB4 vessels (mandatory)
  3. Iba1 microglial cells (optional)

With each image can be provided a *.roi* or *.zip* file containing one or multiple ROI(s):
- if no ROI provided, whole image is analyzed
- if ROI *roiName* provided alone, analysis is performed in ROI not translated along slices
- if ROIs *roiName* and *roiName_2* provided, analysis is performed in ROI translated from *roiName* to *roiName_2* position along slices

### Plugin description

* Detect vessels with Quantile Based Normalization + Cellpose + closing + median filtering
* Compute vessels volume, length, diameter and number of branches, junctions, endpoints, etc
* If microglia channel provided,
  * Detect cells with median filtering + thresholding + closing
  * Compute distance and volume of contact between each microglial cell and its nearest vessel


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *vessels* or *vessels2* models 

### Version history

Version 1 released on October 27, 2023.
