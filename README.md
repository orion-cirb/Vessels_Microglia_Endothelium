# Vessels_Microglia

* **Developed for:** Nicolas
* **Team:** Garel
* **Date:** September 2023
* **Software:** Fiji

### Images description

3D images taken with a x60 objective

3 channels:
  1. IB4 vessels (mandatory)
  2. Lectine vessels (optional)
  3. Iba1 microglial cells (optional)

With each image can be provided a *.roi* or *.zip* file containing one or multiple ROI(s).

### Plugin description

* If 2 vessels channels are provided, add the 2 images
* Detect vessels with Quantile Based Normalization + Cellpose + closing + median filtering
* If microglia channel provided, detect cells with Quantile Based Normalization + Cellpose or median filtering + Otsu thresholding + closing
* Compute distance and volume of contact between each microglial cell and its nearest vessel
* Compute vessels volume, length, diameter and number of branches, junctions, endpoints, etc...

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *vessels*/*vessels2* and *cyto2* models 

### Version history

Version 1 released on September 11, 2023.
