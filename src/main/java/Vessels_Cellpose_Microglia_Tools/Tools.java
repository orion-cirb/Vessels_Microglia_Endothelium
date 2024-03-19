package Vessels_Cellpose_Microglia_Tools;

import Vessels_Cellpose_Microglia_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Vessels_Cellpose_Microglia_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiEnlarger;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import ij.util.ThreadUtil;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Point3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.Measure2Colocalisation;
import mcib3d.geom2.measurements.Measure2Distance;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clijx.bonej.BoneJSkeletonize3D;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.SkeletonResult;


/**
 * @authors Héloïse Monnet & Philippe Mailly
 */
public class Tools {
    
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/Vessels_Cellpose_Microglia";
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    String[] chNames = {"Vessels: ", "Microglia (optional): "};
        
    // Vessels detection
    public String[] vesselSegMethods = {"Cellpose", "Thresholding"};
    public String vesselSegMethod;
    public double minVesselVol = 70;
    // Cellpose method
    public final String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    public final String cellposeModelsPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\.cellpose\\models\\" : System.getProperty("user.home")+"/.cellpose/models/";
    public double cellposeStitchTh = 1;
    public String cellposeModelVessel;
    public int cellposeDiamVessel = 35;
    // Thresholding method
    public String vesselThMethod = "Otsu";
    
    // Microglia segmentation
    public String microThMethod = "Li";
    public double minMicroVol = 20;
    private double roiDilation = 50; // um
    private int microLabel;
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom2.Object3DInt");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        return true;
    }
    
        
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + f);
        }
        Collections.sort(images);
        return(images);
    }
       
    
    /**
     * Find channels name and add None to end of list
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        reader.setId(imageName);
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt = FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
      
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String imagesDir, String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
      
        gd.addMessage("Channels", new Font("Monospace", Font.BOLD, 12), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[n]);
        }
        
        gd.addMessage("Vessels segmentation", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addChoice("Segmentation method", vesselSegMethods, vesselSegMethods[0]);
        gd.addNumericField("Min vessel volume (µm3): ", minVesselVol, 2);
        gd.addMessage("Cellpose", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        String[] models = findCellposeModels();
        gd.addChoice("Cellpose model: ", models, models[0]);
        gd.addMessage("Thresholding", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        String[] thMethods = AutoThresholder.getMethods();
        gd.addChoice("Threshold method: ", thMethods, vesselThMethod); 
        
        gd.addMessage("Microglia segmentation", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addChoice("Threshold method: ", thMethods, microThMethod);        
        gd.addNumericField("Min cell volume (µm3): ", minMicroVol, 2);
        gd.addNumericField("Max cell-vessel distance (µm):", roiDilation, 0);
        
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chChoices = new String[chNames.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();

        vesselSegMethod = gd.getNextChoice();
        minVesselVol = gd.getNextNumber();
        cellposeModelVessel = gd.getNextChoice();
        vesselThMethod = gd.getNextChoice();
        
        microThMethod = gd.getNextChoice();
        minMicroVol = gd.getNextNumber();
        roiDilation = gd.getNextNumber();
        
        if (gd.wasCanceled())
            chChoices = null;
        return(chChoices);
    }
    
        
    /**
     * Find vessel cellpose models in cellpose path
     */
    public String[] findCellposeModels() {
        String[] files = new File(cellposeModelsPath).list();
        ArrayList<String> models = new ArrayList();
        for (String f : files) {
            if (f.contains("vessel")) 
                models.add(f);
        }
        return (models.toArray(new String[0]));
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Preprocess files before sending them to QuantileBasedNormalization plugin
     * @throws Exception
     */
    public void preprocessFiles(ImageProcessorReader reader, ArrayList<String> imageFiles, String processDir, String[] channelNames, String[] channels) throws Exception {
        try {
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setQuiet(true);
                
                // Open and save vessels channel
                int indexCh = ArrayUtils.indexOf(channelNames, channels[0]);
                ImagePlus imgVessels = BF.openImagePlus(options)[indexCh];
                IJ.saveAs(imgVessels, "Tiff", processDir+rootName+"-Vessels.tif");
                closeImage(imgVessels);
            }
        } catch (Exception e) {
            throw e; 
        }
    }
    
        
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        Calibration cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     * Load ROIs, if any provided
     */
    public List<Roi> loadRois(String roiName, ImagePlus img, String imgName) {
        List<Roi> rois = new ArrayList<>();
        
        roiName = new File(roiName+".zip").exists() ? roiName+".zip" : roiName+".roi";
        if (new File(roiName).exists()) {
            RoiManager rm = new RoiManager(false);
            rm.runCommand("Open", roiName);
            List<Roi> roisTemp = Arrays.asList(rm.getRoisAsArray());
            
            for(Roi roi1: roisTemp) {
                if(!roi1.getName().contains("_2")) {
                    boolean translate = false;
                    for(Roi roi2: roisTemp) {
                        if(roi2.getName().equals(roi1.getName()+"_2")) {
                            double x = roi2.getContourCentroid()[0] - roi1.getContourCentroid()[0]; //in pixels
                            double y = roi2.getContourCentroid()[1] - roi1.getContourCentroid()[1]; //in pixels
                            roi1.setProperty("translation", x+"_"+y);
                            rois.add(roi1);
                            translate = true;
                            break;
                        }
                    }
                    if(!translate) {
                        System.out.println("WARNING: ROI " + roi1.getName() + " in image " + imgName + " does not have any associated ROI, it won't be translated across stack");
                        roi1.setProperty("translation", "0_0");
                        rois.add(roi1);
                    }
                }
            }
        } else {
            Roi roi = new Roi(0, 0, img.getWidth() , img.getHeight());
            roi.setName("whole image");
            roi.setProperty("translation", "0_0");
            rois.add(roi);
            System.out.println("WARNING: No ROI file found for image " + imgName + ", entire image is analyzed");
        }

        return(rois);
    }
    
    
    /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     * @throws java.io.IOException
     */
    public ImagePlus vesselSegmentation(ImagePlus img, Calibration cal) throws IOException{
        ImagePlus imgBin = null;
        if(vesselSegMethod == "Cellpose") {
            // Define CellPose settings
            CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModelsPath+cellposeModelVessel, 1, cellposeDiamVessel, cellposeEnvDir);
            settings.setStitchThreshold(cellposeStitchTh);
            settings.useGpu(true);

            // Run CellPose
            ImagePlus imgIn = new Duplicator().run(img);
            CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
            imgBin = cellpose.run();
            
            closeImage(imgIn);
            
        } else if(vesselSegMethod == "Thresholding") {
            ImagePlus imgMed = medianFilter(img, true, 2, 0);
            ImagePlus imgDOG = DOG(imgMed, 5, 10);
            imgBin = threshold(imgDOG, vesselThMethod);

            closeImage(imgMed);
            closeImage(imgDOG);
        }
        
        // Remove small objects and connect remaining ones
        ImagePlus imgClose = closingFilter(imgBin, 6, 1);
        ImagePlus imgOut = medianFilter(imgClose, true, 1, 0);
        imgOut.setCalibration(cal);

        // Get population of detections
        Objects3DIntPopulation pop = getPopFromImage(imgOut);
        System.out.println("Nb objects detected: " + pop.getNbObjects());
        popFilterSize(pop, minVesselVol, Double.MAX_VALUE);
        System.out.println("Nb objects remaining after size filtering: " + pop.getNbObjects());
        
        // Draw population in image
        ImageHandler imhMask = ImageHandler.wrap(img).createSameDimensions();
        for (Object3DInt obj: pop.getObjects3DInt())
            obj.drawObject(imhMask, 255);

        closeImage(imgBin);
        closeImage(imgClose);
        closeImage(imgOut);
        return(imhMask.getImagePlus());
    }
    
    
    /**
     * Closing filtering using CLIJ2
     */ 
    private ImagePlus closingFilter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DBox(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);

       ClearCLBuffer imgCLMin = clij2.create(imgCLMax);
       clij2.minimum3DBox(imgCLMax, imgCLMin, sizeXY, sizeXY, sizeZ);
       ImagePlus imgMin = clij2.pull(imgCLMin);
       
       clij2.release(imgCL);
       clij2.release(imgCLMax);
       clij2.release(imgCLMin);
       return(imgMin);
    }
    
    
    /**
     * 3D median filtering using CLIJ2
     */ 
    public ImagePlus medianFilter(ImagePlus img, boolean sliceBySlice, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCL = clij2.push(img); 
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        if (sliceBySlice)
            clij2.median3DSliceBySliceSphere(imgCL, imgCLMed, sizeXY, sizeXY);
        else
            clij2.median3DSphere(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        ImagePlus imgMed = clij2.pull(imgCLMed);
        clij2.release(imgCL);
        clij2.release(imgCLMed);
        return(imgMed);
    }
    
            
    /**
     * Difference of Gaussians using CLIJ
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, size1, size1, size2, size2);
        clij2.release(imgCL);
        ImagePlus imgDOG = clij2.pull(imgCLDOG); 
        clij2.release(imgCLDOG);
        return(imgDOG);
    }

    
    /**
     * Automatic thresholding using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
   
    /**
     * Return population of 3D objects population from binary image
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        labels.closeImagePlus();
        return(pop);
    }
      
    
    /**
     * Remove object with size < min and size > max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
      
    
    /**
     * Segment microglia with median filtering + Otsu thresholding + closing filtering
     */
    public Objects3DIntPopulation cellsSegmentation(ImagePlus imgIn, Calibration cal) {
        ImagePlus imgMed = medianFilter(imgIn, false, 2, 1);
        ImagePlus imgBin = threshold(imgMed, microThMethod);
        ImagePlus imgClose = closingFilter(imgBin, 2, 2);
        imgClose.setCalibration(cal);
        
        Objects3DIntPopulation pop = getPopFromImage(imgClose);
        System.out.println("Nb objects detected: " + pop.getNbObjects());
        popFilterSize(pop, minMicroVol, Double.MAX_VALUE);
        System.out.println("Nb objects remaining after size filtering: " + pop.getNbObjects());
        
        closeImage(imgMed);
        closeImage(imgBin);
        closeImage(imgClose);
        return(pop);
    }
    
    
    /** 
     * Compute parameters and save results for each ROI
     * @throws java.io.IOException
     */
    public void saveResults(ImagePlus imgVesselMask, Objects3DIntPopulation microPop, List<Roi> rois,
            ImagePlus imgVessels, ImagePlus imgMicro, Calibration cal, BufferedWriter microResults, BufferedWriter globalResults, 
            String rootName, String outDir) throws IOException {
        
        ImageHandler imhVessels = ImageHandler.wrap(imgVessels).createSameDimensions();
        ImageHandler imhSkel = imhVessels.createSameDimensions();
        ImageHandler imhMicro = imhVessels.createSameDimensions();
        microLabel = 1;
        
        for (Roi roi: rois) {
            print("- Computing parameters and saving results for ROI " + roi.getName() + " -");
            
            // Write ROI parameters in global results file
            boolean roiTranslation = roi.getProperty("translation").equals("0_0")? false: true;
            imgVessels.setRoi(roi);
            double roiVol = imgVessels.getStatistics().area * imgVessels.getNSlices() * cal.pixelDepth;
            imgVessels.deleteRoi();
            globalResults.write(rootName+"\t"+cal.pixelWidth+"\t"+roi.getName()+"\t"+roiTranslation+"\t"+roiVol);
            globalResults.flush();

            // Get vessels in translated ROI
            ImagePlus imgVesselMaskRoi = clearOutsideRoi(imgVesselMask, roi, false, cal);
            double vesselVol = new MeasureVolume(new Object3DInt(ImageHandler.wrap(imgVesselMaskRoi))).getVolumeUnit();
            double vesselDensity = vesselVol / roiVol * 1e6;

            // Compute vessels skeleton
            ImagePlus imgVesselSkel = skeletonize3D(imgVesselMaskRoi, cal);
            Object3DInt objVesselSkel = new Object3DInt(ImageHandler.wrap(imgVesselSkel));
            IJ.run(imgVesselSkel, "8-bit","");
            
            // Get vessels in translated dilated ROI
            ImagePlus imgVesselMaskRoiDil = clearOutsideRoi(imgVesselMask, roi, true, cal);
            Object3DInt objVesselMaskRoiDil = new Object3DInt(ImageHandler.wrap(imgVesselMaskRoiDil));
            
            // Compute vessels distance map
            ImageFloat vesselDistMap = localThickness3D(imgVesselMaskRoiDil, false, cal);

            // Write vessels parameters in global results file
            writeVesselParams(imgVesselSkel, vesselDistMap.getImagePlus(), vesselVol, vesselDensity, globalResults);

            Objects3DIntPopulation microPopInRoi = new Objects3DIntPopulation();
            if (imgMicro != null) {
                // Get microglia in translated ROI
                microPopInRoi = getPopInsideRoi(microPop, imgMicro, roi);
                
                // Compute vessels inverse distance maps
                ImageFloat vesselDistMapInv = localThickness3D(imgVesselMaskRoiDil, true, cal);

                // Write microglia parameters in micro and global results files
                writeMicroParams(microPopInRoi, objVesselMaskRoiDil, objVesselSkel, vesselDistMap, vesselDistMapInv, vesselVol, cal, rootName, roi.getName(), microResults, globalResults);
            
                vesselDistMapInv.closeImagePlus();
            }
            globalResults.write("\n");
            globalResults.flush();

            // Draw results
            System.out.println("Drawing results...");
            objVesselMaskRoiDil.drawObject(imhVessels, 128);
            objVesselSkel.drawObject(imhSkel, 255);
            if (microPopInRoi.getNbObjects() != 0)
                microPopInRoi.drawInImage(imhMicro);

            closeImage(imgVesselMaskRoi);
            closeImage(imgVesselSkel);
            closeImage(imgVesselMaskRoiDil);
            vesselDistMap.closeImagePlus();
        }
        
        drawResults(imhVessels, imhSkel, imhMicro, imgVessels, imgMicro, cal, rootName, outDir);
        imhVessels.closeImagePlus();
        imhSkel.closeImagePlus();
        imhMicro.closeImagePlus();
    }  
    
    
    /**
     * Clear image outside ROI translated across stack 
     */
    private ImagePlus clearOutsideRoi(ImagePlus imgIn, Roi roi, boolean dilate, Calibration cal) {
        ImagePlus img = new Duplicator().run(imgIn);
                
        int nbSlices = img.getNSlices();
        String[] translation = roi.getProperty("translation").split("_");
        double x = Double.valueOf(translation[0]);
        double y = Double.valueOf(translation[1]);
        double x_step = x / (nbSlices-1);
        double y_step = y / (nbSlices-1);
        
        img.setRoi(dilate ? RoiEnlarger.enlarge(roi, roiDilation/cal.pixelWidth) : roi);
        for(int i=1; i <= nbSlices; i++) {
            img.setPosition(i);
            if(i > 1)
                IJ.run(img, "Translate... ", "x="+x_step+" y="+y_step); // in pixels
            IJ.run(img, "Clear Outside", "slice");
        }
        IJ.run(img, "Translate... ", "x="+(-1*x)+" y="+(-1*y)); // in pixels
        img.deleteRoi();
        
        return(img);
    }
    
    
    /**
     * Compute (inverse) distance map
     */
    public ImageFloat localThickness3D(ImagePlus img, boolean inverse, Calibration cal) {
        System.out.println("Computing vessels (inverted) distance map...");
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0, inverse, ThreadUtil.getNbCpus());
        edt.setCalibration(cal);
        return(edt);
    }
    
    
    /**
     * Skeletonize 3D with CLIJ2
     */
    public ImagePlus skeletonize3D(ImagePlus img, Calibration cal) {
        System.out.println("Computing vessels skeleton...");
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLSkel = clij2.create(imgCL);
        new BoneJSkeletonize3D().bonejSkeletonize3D(clij2,imgCL, imgCLSkel);
        ImagePlus imgSkel = clij2.pull(imgCLSkel);
        clij2.release(imgCL);
        clij2.release(imgCLSkel);
        
        imgSkel.setCalibration(cal);
        return(imgSkel);
    }
    
    
    /**
     * Compute skeleton parameters and write all vessels parameters in results file (length, number of branches, junctions, endpoints...)
     */ 
    private void writeVesselParams(ImagePlus skel, ImagePlus distMap, double vesselVol, double vesselDensity, 
            BufferedWriter globalResults) throws IOException {       
        System.out.println("Computing vessels parameters...");
        if(vesselVol == 0) {
            globalResults.write("\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");
            globalResults.flush();
        } else {
            AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
            analyzeSkeleton.setup("", skel);
            // Params: int pruneIndex, boolean pruneEnds, boolean shortPath, ImagePlus origIP, boolean silent, boolean verbose
            SkeletonResult skeletonResult = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, false, null, true, false);

            double[] branchLengths = skeletonResult.getAverageBranchLength();
            int[] branchNumbers = skeletonResult.getBranches();
            int totalLength = 0;
            for (int i = 0; i < branchNumbers.length; i++)
                totalLength += branchNumbers[i] * branchLengths[i];
            double meanLength = StatUtils.mean(branchLengths);
            double lengthLongestBranch = StatUtils.max(skeletonResult.getMaximumBranchLength());
            int nbBranches = IntStream.of(branchNumbers).sum();
            int nbJunctions = IntStream.of(skeletonResult.getJunctions()).sum();
            int nbEndpoints = IntStream.of(skeletonResult.getEndPoints()).sum();

            ArrayList<Point> voxels = skeletonResult.getListOfSlabVoxels();
            voxels.addAll(skeletonResult.getListOfJunctionVoxels());
            voxels.addAll(skeletonResult.getListOfEndPoints());
            DescriptiveStatistics diameters = new DescriptiveStatistics();
            for (Point pt: voxels) {
                distMap.setZ(pt.z);
                double diam = distMap.getProcessor().getPixelValue(pt.x, pt.y)*2;
                if (diam > 0) diameters.addValue(diam);
            }

            globalResults.write("\t"+vesselVol+"\t"+vesselDensity+"\t"+totalLength+"\t"+meanLength+"\t"+
                    lengthLongestBranch+"\t"+nbBranches+"\t"+nbJunctions+"\t"+nbEndpoints+"\t"+
                    diameters.getMean()+"\t"+diameters.getStandardDeviation()+"\t"+
                    diameters.getMin()+"\t"+diameters.getMax());
            globalResults.flush();
        }
    }
    
    
    /**
     * Get population of objects with their centroid into ROI translated across stack 
     */
     private Objects3DIntPopulation getPopInsideRoi(Objects3DIntPopulation popIn, ImagePlus img, Roi roi) {
        Objects3DIntPopulation popOut = new Objects3DIntPopulation();
        
        int nbSlices = img.getNSlices();
        String[] translation = roi.getProperty("translation").split("_");
        double x = Double.valueOf(translation[0]);
        double y = Double.valueOf(translation[1]);
        double x_step = x / (nbSlices-1);
        double y_step = y / (nbSlices-1);
        
        for(Object3DInt obj: popIn.getObjects3DInt()) {
            Point3D centroid = new MeasureCentroid(obj).getCentroidAsPoint();
            if(centroid.getRoundZ() > 1)    
                roi.translate(centroid.getRoundZ()*x_step, centroid.getRoundZ()*y_step);
            if(roi.contains(centroid.getRoundX(), centroid.getRoundY()))
               popOut.addObject(obj);
            if(centroid.getRoundZ() > 1) 
                roi.translate(-1*centroid.getRoundZ()*x_step, -1*centroid.getRoundZ()*y_step);
        }
        return(popOut);
    }
     
     
    /**
     * Compute and write microglia parameters in results files
     */ 
    private void writeMicroParams(Objects3DIntPopulation microPop, Object3DInt vesselObj, Object3DInt vesselSkelObj, 
            ImageFloat vesselDistMap, ImageFloat vesselDistMapInv, double vesselVol, Calibration cal, 
            String rootName, String roiName, BufferedWriter microResults, BufferedWriter globalResults) throws IOException {
        System.out.println("Computing microglia parameters...");
        int nbVAM=0, nbVTM=0, nbVDM = 0;
        for (Object3DInt micro: microPop.getObjects3DInt()) {
            double microVol = new MeasureVolume(micro).getVolumeUnit();
            microResults.write(rootName+"\t"+roiName+"\t"+microLabel+"\t"+microVol);
            micro.setLabel(microLabel);
            microLabel++;
            if(vesselVol == 0) {
                microResults.write("\t"+0+"\t"+Double.NaN+"\t"+Double.NaN+"\t"+Double.NaN+"\n");

                nbVDM++;
            } else {
                double borderDist = new Measure2Distance(micro, vesselObj).getValue(Measure2Distance.DIST_BB_UNIT);
                double vesselDiam = 2*vesselDistMap.getPixel(new Measure2Distance(micro, vesselSkelObj).getBorder2Pix());
                double colocVol = getColocVol(micro, vesselObj, cal);
                double centroidDist = vesselDistMapInv.getPixel(new MeasureCentroid​(micro).getCentroidAsPoint());
                microResults.write("\t"+colocVol+"\t"+centroidDist+"\t"+borderDist+"\t"+vesselDiam+"\n");

                if(colocVol == 0) nbVDM++;
                else if(colocVol != 0 && centroidDist == 0) nbVAM++;
                else if(colocVol != 0 && centroidDist != 0) nbVTM++;
            }
            microResults.flush();
        }
        globalResults.write("\t"+microPop.getNbObjects()+"\t"+nbVAM+"\t"+nbVTM+"\t"+nbVDM);
        globalResults.flush();
    }
     
    
    /**
     * Get colocalization volume between microglia and vessel 
     */
    private double getColocVol(Object3DInt micro, Object3DInt vessel, Calibration cal) {
        Measure2Colocalisation coloc = new Measure2Colocalisation(micro, vessel);
        double colocVol = coloc.getValue(Measure2Colocalisation.COLOC_VOLUME);
        double pixVol = cal.pixelHeight*cal.pixelWidth*cal.pixelDepth;
        return(colocVol*pixVol);
    }
    
    
    public void drawResults(ImageHandler imhVessels, ImageHandler imhSkel, ImageHandler imhMicro, ImagePlus imgVessels, ImagePlus imgMicro, Calibration cal, String rootName, String outDir)  {
        IJ.run(imhVessels.getImagePlus(), "Red", "");
        IJ.run(imhSkel.getImagePlus(), "Blue", "");
        IJ.run(imhMicro.getImagePlus(), "Green", "");
        
        ImagePlus[] imgColors = {imhVessels.getImagePlus(), null, imhSkel.getImagePlus(), imgVessels, null};
        if (imgMicro != null) {
            imgColors[1] = imhMicro.getImagePlus();
            imgColors[4] = imgMicro;
        }
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(cal);
        FileSaver imgObjectsFile = new FileSaver(imgObjects);
        imgObjectsFile.saveAsTiff(outDir+rootName+".tif");
        closeImage(imgObjects);
    }
   
}
