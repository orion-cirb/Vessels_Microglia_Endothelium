package Vessels_Microglia_Endothelium_Tools;

import Vessels_Microglia_Endothelium_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Vessels_Microglia_Endothelium_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiEnlarger;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.util.ThreadUtil;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
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
import mcib3d.geom2.tracking.TrackingAssociation;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clijx.bonej.BoneJSkeletonize3D;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.SkeletonResult;


/**
 * @authors Héloïse Monnet & Philippe Mailly
 */
public class Tools {
    
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/Vessels_Microglia_Endothelium";
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    private String[] chNames = {"Vessels: ", "Microglia (optional): ", "Endothelial nuclei (optional): "};
        
    // Vessels detection
    public boolean vesselNormalization = true;
    private String[] vesselSegMethods = {"Thresholding", "Cellpose"};
    public String vesselSegMethod;
    private double minVesselVol = 70; // um3
    private double minVesselLength = 10; // um
    // Cellpose method
    private final String cellposeEnvPath = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    private final String cellposeModelsPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\.cellpose\\models\\" : System.getProperty("user.home")+"/.cellpose/models/";
    public String cellposeModelVessel;
    private int cellposeDiamVessel = 35; // pix
    private double cellposeStitchThVessel = 1;
    // Thresholding method
    public String vesselThMethod = "RenyiEntropy";
    
    // Microglia segmentation
    public String microThMethod = "Li";
    private double microMinVol = 20; // um3
    private double roiDilation = 50; // um
    
    // Endothelial nuclei detection
    private double endoMinVol = 20; // um3
    // Omnipose
    private String omniposeEnvPath = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"omnipose" : "/opt/miniconda3/envs/omnipose";
    private String omniposeModelsPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\" : System.getProperty("user.home")+"/.cellpose/models/";
    private String omniposeModelEndo = "cyto2_omnitorch_0";
    private int omniposeDiamEndo = 25;
    // Masks stitching
    private float maxLabel;
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
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
     * Get images with given extension in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
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
     * Get channels name and add None to the end of the list
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
     * Display dialog box
     */
    public String[] dialog(String imagesDir, String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 70, 0);
        gd.addImage(icon);
      
        gd.addMessage("Channels", new Font("Monospace", Font.BOLD, 12), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[n]);
        }
        
        gd.addMessage("Vessels segmentation", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addNumericField("Min vessel volume (µm3): ", minVesselVol, 2);
        gd.addNumericField("Min branch length (µm): ", minVesselLength, 2);
        gd.addCheckbox("Quantile based normalization", vesselNormalization);
        gd.addChoice("Segmentation method", vesselSegMethods, vesselSegMethods[0]);
        gd.addMessage("Thresholding", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        String[] thMethods = AutoThresholder.getMethods();
        gd.addChoice("Threshold method: ", thMethods, vesselThMethod);
        gd.addMessage("Cellpose", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        String[] models = findCellposeModels();
        gd.addChoice("Cellpose model: ", models, models[0]);
        
        gd.addMessage("Microglia segmentation", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addChoice("Threshold method: ", thMethods, microThMethod);        
        gd.addNumericField("Min cell volume (µm3): ", microMinVol, 2);
        gd.addNumericField("ROI dilation (µm):", roiDilation, 0);
        
        gd.addMessage("Endothelial nuclei segmentation", new Font("Monospace", Font.BOLD, 12), Color.blue); 
        gd.addNumericField("Min nucleus volume (µm3): ", endoMinVol, 2);
        
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chChoices = new String[chNames.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();

        minVesselVol = gd.getNextNumber();
        minVesselLength = gd.getNextNumber();
        vesselNormalization = gd.getNextBoolean();
        vesselSegMethod = gd.getNextChoice();
        vesselThMethod = gd.getNextChoice();
        cellposeModelVessel = gd.getNextChoice();
        
        microThMethod = gd.getNextChoice();
        microMinVol = gd.getNextNumber();
        roiDilation = gd.getNextNumber();
        
        endoMinVol = gd.getNextNumber();
        
        if (gd.wasCanceled())
            chChoices = null;
        return(chChoices);
    }
    
        
    /**
     * Get vessel Cellpose models in Cellpose models directory
     */
    private String[] findCellposeModels() {
        String[] files = new File(cellposeModelsPath).list();
        ArrayList<String> models = new ArrayList();
        for (String f : files) {
            if (f.contains("vessel")) 
                models.add(f);
        }
        return (models.toArray(new String[0]));
    }
    
    
    /**
     * Write headers in results files
     */
    public void writeHeaders(String[] channels, BufferedWriter globalResults, BufferedWriter vesselResults, BufferedWriter microResults) throws IOException {         
        globalResults.write("Image name\tImage XY calibration (µm)\tROI name\tROI translation\tROI volume (µm3)"
                + "\tVessels total volume (µm3)\tVessels density\tVessels total length (µm)"
                + "\tNb branches\tNb junctions\tVessels mean diameter (µm)\tVessels std diameter (µm)");
            if(channels[1] != "None") 
                globalResults.write("\tNb microglia\tNb VAM\tNb VTM\tNb VDM");
            if(channels[2] != "None") 
                globalResults.write("\tNb endothelial nuclei\tEndothelial nuclei density\tNb endothelial nuclei/Vessels total length");
        globalResults.write("\n");
        globalResults.flush();

        vesselResults.write("Image name\tROI name\tBranch length (µm)\tMean diameter (µm)\tDiameter std (µm)\tMin diameter (µm)\tMax diameter (µm)"
                + "\tV1 x (µm)\tV1 y (µm)\tV1 z (µm)\tV2 x (µm)\tV2 y (µm)\tV2 z (µm)\n");
        vesselResults.flush();
        
        if(channels[1] != "None")  {
            microResults.write("Image name\tROI name\tCell ID\tCell volume (µm3)\tCell volume coloc with vessel (µm3)"
                    + "\tCell centroid distance to closest vessel (µm)\tCell border distance to closest vessel (µm)"
                    + "\tClosest vessel diameter (µm)\n");
            microResults.flush();
        }
    }
    
    
    /**
     * Save images specific channel before sending it to QuantileBasedNormalization plugin
     * @throws Exception
     */
    public void saveChannel(String dir, ArrayList<String> imageFiles, int channelIndex, String extension) throws Exception {
        try {
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setQuiet(true);
                
                // Open and save vessels channel
                ImagePlus imgVessels = BF.openImagePlus(options)[channelIndex];
                IJ.saveAs(imgVessels, "Tiff", dir+rootName+extension);
                closeImage(imgVessels);
            }
        } catch (Exception e) {
            throw e; 
        }
    }
    
    
    /**
     * Delete images specific channel after it was sent to QuantileBasedNormalization plugin
     * @throws Exception
     */
    public void deleteChannel(String dir, ArrayList<String> imageFiles, String extension) throws Exception {
        try {
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                new File(dir+rootName+extension).delete();
            }
        } catch (Exception e) {
            throw e; 
        }
    }
    
        
    /**
     * Get image calibration
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
     * Detect 3D vessels in a Z-stack with 2 different methods:
     * - Cellpose 2D applied slice by slice 
     * - Median filtering + DoG filtering + automatic thresholding applied slice by slice 
     * @throws java.io.IOException
     */
    public ImagePlus vesselSegmentation(ImagePlus img, Calibration cal) throws IOException{
        ImagePlus imgBin = null;
        if(vesselSegMethod == "Cellpose") {
            // Define CellPose settings
            CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModelsPath+cellposeModelVessel, 1, cellposeDiamVessel, cellposeEnvPath);
            settings.setStitchThreshold(cellposeStitchThVessel);
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
        ImagePlus imgOut = medianFilter(imgClose, false, 2, 1);
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
     * Segment microglia with median filtering + thresholding + closing filtering
     */
    public Objects3DIntPopulation microSegmentation(ImagePlus img, Calibration cal) {
        ImagePlus imgMed = medianFilter(img, false, 2, 1);
        ImagePlus imgBin = threshold(imgMed, microThMethod);
        ImagePlus imgClose = closingFilter(imgBin, 2, 2);
        imgClose.setCalibration(cal);
        
        // Get population of objects
        Objects3DIntPopulation pop = getPopFromImage(imgClose);
        System.out.println("Nb objects detected: " + pop.getNbObjects());
        
        // Filter objects
        popFilterSize(pop, microMinVol, Double.MAX_VALUE);
        System.out.println("Nb objects remaining after size filtering: " + pop.getNbObjects());
        
        closeImage(imgMed);
        closeImage(imgBin);
        closeImage(imgClose);
        return(pop);
    }
    
    
    /**
     * Segment endothelial nuclei with median filtering + Omnipose
     */
    public Objects3DIntPopulation endoSegmentation(ImagePlus img, Calibration cal) {
        int nSlices = img.getDimensions()[3];
       
        // Median filter
        ImagePlus imgMed = medianFilter(img, false, 1, 1);
        imgMed.setDimensions​(1, 1, nSlices);
        
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModelsPath+omniposeModelEndo, 1, omniposeDiamEndo, omniposeEnvPath);
        settings.setVersion("0.7");
        settings.setCluster(true);
        settings.setOmni(true);
        settings.useGpu(true);

        // Run Omnipose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgMed);
        PrintStream console = System.out;
        System.setOut(new NullPrintStream());
        ImagePlus imgBin = cellpose.run();
        System.setOut(console);
        imgBin.setDimensions​(1, nSlices, 1);
        imgBin.setCalibration(cal);
        
        // Get population of objects
        ImagePlus imgStitch = stitch3D(imgBin);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgStitch));
        System.out.println("Nb objects detected: " + pop.getNbObjects());

        // Filter objects
        popFilterOneZ(pop);
        popFilterSize(pop, endoMinVol, Double.MAX_VALUE);
        System.out.println("Nb objects remaining after filtering: " + pop.getNbObjects());

        closeImage(imgMed);
        closeImage(imgBin);
        closeImage(imgStitch);
        return(pop);
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
    private ImagePlus medianFilter(ImagePlus img, boolean sliceBySlice, double sizeXY, double sizeZ) {
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
    private ImagePlus DOG(ImagePlus img, double size1, double size2) {
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
    private ImagePlus threshold(ImagePlus img, String thMed) {
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
    private Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageInt labels = new ImageLabeller().getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        labels.closeImagePlus();
        return(pop);
    }
    
    
    /**
     * Stitch 2D masks into 3D volumes in a stack
     */
    private ImagePlus stitch3D(ImagePlus img) {
        maxLabel = 0;
        IJ.run(img, "Select None", "");
        
        ImagePlus[] slices = new ImagePlus[img.getNSlices()];
        slices[0] = img.crop(1+"-"+1);
        for (int i = 1; i < img.getNSlices(); i++) {
            ImagePlus nextSlice = img.crop((i+1)+"-"+(i+1));
            slices[i] = associate(slices[i-1], nextSlice);
            closeImage(nextSlice);            
        }
        
        ImagePlus imgStitch = new Concatenator().concatenate(slices, false);
        imgStitch.setDimensions(1, imgStitch.getNFrames(), 1);
        return(imgStitch);
    }
    
    
    /** 
     * Associate 2D masks of slice z-1 with 2D masks of slice z
     */
    private ImagePlus associate(ImagePlus imp1, ImagePlus imp2) {
        TrackingAssociation association = new TrackingAssociation(ImageInt.wrap(imp1), ImageInt.wrap(imp2), 0, 0.1);
        association.setMaxLabel(maxLabel);
        ImageHandler imp2Associated = association.getTrackedImage();
        maxLabel = association.getMaxLabel();
        return(imp2Associated.getImagePlus());
    }
      
    
    /**
     * Remove object with size < min and size > max
     */
    private void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    /**
     * Remove objects in population with only one plan
     */
    private void popFilterOneZ(Objects3DIntPopulation pop) {
        pop.getObjects3DInt().removeIf(p -> (p.getObject3DPlanes().size() == 1));
        pop.resetLabels();
    }
  
    
    /**
     * Skeletonize 3D with CLIJ2
     */
    public ImagePlus skeletonize3D(ImagePlus img, Calibration cal) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLSkel = clij2.create(imgCL);
        new BoneJSkeletonize3D().bonejSkeletonize3D(clij2,imgCL, imgCLSkel);
        ImagePlus imgSkel = clij2.pull(imgCLSkel);
        clij2.release(imgCL);
        clij2.release(imgCLSkel);
        
        IJ.run(imgSkel, "8-bit","");
        imgSkel.setCalibration(cal);
        return(imgSkel);
    }
    
    
    /**
     * Prune skeleton branches with length smaller than threshold
     * https://imagej.net/plugins/analyze-skeleton/
     */
    public ImagePlus pruneSkeleton(ImagePlus image) {
        // Analyze skeleton
        AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();
        skel.setup("", image);
        SkeletonResult skelResult = skel.run(AnalyzeSkeleton_.NONE, false, false, null, true, false);

        // Create copy of input image
        ImagePlus prunedImage = image.duplicate();
        if(skelResult.getBranches() != null) {
            ImageStack outStack = prunedImage.getStack();

            // Get graphs (one per skeleton in the image)
            Graph[] graphs = skelResult.getGraph();

            // Get list of end-points
            ArrayList<Point> endPoints = skelResult.getListOfEndPoints();

            for(Graph graph: graphs) {
                ArrayList<Edge> listEdges = graph.getEdges();

                // Go through all branches and remove branches under threshold in duplicate image
                for(Edge e: listEdges) {
                    ArrayList<Point> p1 = e.getV1().getPoints();
                    boolean v1End = endPoints.contains(p1.get(0));
                    ArrayList<Point> p2 = e.getV2().getPoints();
                    boolean v2End = endPoints.contains(p2.get(0));
                    // If any of the vertices is end-point 
                    if(v1End || v2End) {
                        if(e.getLength() < minVesselLength) { // in microns
                            if(v1End)
                                outStack.setVoxel(p1.get(0).x, p1.get(0).y, p1.get(0).z, 0);
                            if(v2End)
                                outStack.setVoxel(p2.get(0).x, p2.get(0).y, p2.get(0).z, 0);
                            for(Point p: e.getSlabs())
                                outStack.setVoxel(p.x, p.y, p.z, 0);
                        }
                    }
                }
            }
        }
        
        return(prunedImage);
    }
    
    
    /** 
     * Compute parameters and save results for each ROI
     * @throws java.io.IOException
     */
    public void saveResults(List<Roi> rois, ImagePlus imgVesselMask, ImagePlus imgVesselSkel, ImageFloat vesselDistMap, ImageFloat vesselDistMapInv, 
                            Objects3DIntPopulation microPop, Objects3DIntPopulation endoPop, ImagePlus imgVessels, ImagePlus imgMicro, 
                            ImagePlus imgEndo, Calibration cal, BufferedWriter vesselResults, BufferedWriter microResults, BufferedWriter globalResults, 
                            String imgName, String dirName) throws IOException {        
        
        ImageHandler imhVessels = ImageHandler.wrap(imgVessels).createSameDimensions();
        ImageHandler imhMicro = imhVessels.createSameDimensions();
        ImageHandler imhEndo = imhVessels.createSameDimensions();
        
        ImageStack stackTagSkel = new ImageStack(imgVessels.getWidth(), imgVessels.getHeight());
        for (int z = 0; z < imgVessels.getNSlices(); z++)
            stackTagSkel.addSlice(new ByteProcessor(imgVessels.getWidth(), imgVessels.getHeight()));
        
        AtomicInteger microLabel = new AtomicInteger(1);

        for (Roi roi: rois) {
            print("- Computing parameters and saving results for ROI " + roi.getName() + " -");
            
            // VESSELS
            // Get vessels in (translated) non-dilated ROI
            ImagePlus imgVesselMaskRoi = clearOutsideRoi(imgVesselMask, roi, false, cal);
            Object3DInt objVesselMaskRoi = new Object3DInt(ImageHandler.wrap(imgVesselMaskRoi));
            double vesselVol = new MeasureVolume(objVesselMaskRoi).getVolumeUnit();
            if (imgMicro == null) objVesselMaskRoi.drawObject(imhVessels, 128);
            closeImage(imgVesselMaskRoi);
            
            // Begin to write parameters in global results file
            double roiVol = computeRoiVolume(roi, imgVessels, cal);
            globalResults.write(imgName+"\t"+cal.pixelWidth+"\t"+roi.getName()+"\t"+!roi.getProperty("translation").equals("0_0")+
                                "\t"+roiVol+"\t"+vesselVol+"\t"+vesselVol/roiVol*1e6);
            globalResults.flush();
            
            // Write vessels skeleton parameters in global results file
            double vesselTotalLength = saveVesselResultsInRoi(roi, imgVesselSkel, vesselDistMap, stackTagSkel, 
                                                              globalResults, vesselResults, cal, imgName);

            // MICROGLIA
            if (imgMicro != null) {
                // Get vessels in (translated) dilated ROI
                ImagePlus imgVesselMaskRoiDil = clearOutsideRoi(imgVesselMask, roi, true, cal);
                new Object3DInt(ImageHandler.wrap(imgVesselMaskRoiDil)).drawObject(imhVessels, 128);
                saveMicroResultsInRoi(roi, microPop, imgMicro, imhMicro, microLabel, imgVesselMaskRoiDil, imgVesselSkel,
                                      vesselDistMap, vesselDistMapInv, cal, vesselVol, microResults, globalResults, imgName, roi.getName());
                closeImage(imgVesselMaskRoiDil);
            }
            
            // ENDOTHELIAL NUCLEI
            if (imgEndo != null)
                saveEndoResultsInRoi(roi, endoPop, imgEndo, imhEndo, roiVol, vesselTotalLength, globalResults);
            
            globalResults.write("\n");
            globalResults.flush();
        }
        
        // Save drawings
        saveImages(stackTagSkel, imhVessels, imhMicro, imhEndo, imgVessels, imgMicro, imgEndo, cal, imgName, dirName);
        
        imhVessels.closeImagePlus();
        imhMicro.closeImagePlus();
        imhEndo.closeImagePlus();
        
    }
    
    
    /**
     * Compute ROI volume 
     */
    private double computeRoiVolume(Roi roi, ImagePlus img, Calibration cal) {
        img.setRoi(roi);
        double roiVolume = img.getStatistics().area * img.getNSlices() * cal.pixelDepth;
        img.deleteRoi();
        return(roiVolume);
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
    public ImageFloat distanceMap3D(ImagePlus img, boolean inverse, Calibration cal) {
        img.setCalibration(cal);
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0, inverse, ThreadUtil.getNbCpus());
        return(edt);
    }
    
    
    /**
     * Compute vessels parameters in ROI and write them in global results file
     */
    private double saveVesselResultsInRoi(Roi roi, ImagePlus imgVesselSkel, ImageFloat distMap, ImageStack stackTagSkel,
                                          BufferedWriter globalResults, BufferedWriter vesselResults, Calibration cal, String imgName) throws IOException {
        System.out.println("Computing vessels parameters...");
        
        // Get skeleton in (translated) non-dilated ROI
        ImagePlus imgVesselSkelRoi = clearOutsideRoi(imgVesselSkel, roi, false, cal);
        AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();
        skel.setup("", imgVesselSkelRoi);
        SkeletonResult skelResult = skel.run(AnalyzeSkeleton_.NONE, false, false, null, true, false);
        closeImage(imgVesselSkelRoi);
        
        // Save parameters in global results file
        double totalLength = 0;
        if(skelResult.getBranches() == null) {
            globalResults.write("\t0\t0\t0\t0\t0");
            globalResults.flush(); 
        } else {
            double[] branchLengths = skelResult.getAverageBranchLength();
            int[] branchNumbers = skelResult.getBranches();
            for (int i = 0; i < branchNumbers.length; i++)
                totalLength += branchNumbers[i] * branchLengths[i];
            int nbBranches = IntStream.of(branchNumbers).sum();
            int nbJunctions = IntStream.of(skelResult.getJunctions()).sum();
            
            DescriptiveStatistics diameters = new DescriptiveStatistics();
            for (Point pt: skelResult.getListOfSlabVoxels())
                diameters.addValue(distMap.getPixel(pt.x, pt.y, pt.z)*2);
                
            globalResults.write("\t"+totalLength+"\t"+nbBranches+"\t"+nbJunctions
                                +"\t"+diameters.getMean()+"\t"+diameters.getStandardDeviation());
            globalResults.flush();
        }
        
        // Save parameters in vessel results file
        for(Graph graph: skelResult.getGraph()) {
            ArrayList<Edge> edges = graph.getEdges();
            for(Edge e: edges) {
                DescriptiveStatistics diams = new DescriptiveStatistics();
                for (Point pt: e.getSlabs())
                    diams.addValue(distMap.getPixel(pt.x, pt.y, pt.z)*2);

                Point v1 = e.getV1().getPoints().get(0);
                Point v2 = e.getV2().getPoints().get(0);
                vesselResults.write(imgName+"\t"+roi.getName()+"\t"+e.getLength()+"\t"+diams.getMean()+
                                    "\t"+diams.getStandardDeviation()+"\t"+diams.getMin()+"\t"+diams.getMax()+
                                    "\t"+v1.x*cal.pixelWidth+"\t"+v1.y*cal.pixelHeight+"\t"+v1.z*cal.pixelDepth+
                                    "\t"+v2.x*cal.pixelWidth+"\t"+v2.y*cal.pixelHeight+"\t"+v2.z*cal.pixelDepth+"\n"); 
                vesselResults.flush();
            }
        }
        
        tagImage(stackTagSkel, skelResult);
        
        return(totalLength);
    }

    
    /**
     * Tag skeleton dividing the voxels between slabs, end points and junctions
     */
    private ImageStack tagImage(ImageStack stack, SkeletonResult skeletonResult) {        
        for(Point p: skeletonResult.getListOfSlabVoxels())
            stack.setVoxel(p.x, p.y, p.z, 40);
        for(Point p: skeletonResult.getListOfEndPoints())
            stack.setVoxel(p.x, p.y, p.z, 40);
        for(Point p: skeletonResult.getListOfJunctionVoxels())
            stack.setVoxel(p.x, p.y, p.z, 200);
        
        return(stack);
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
    private void saveMicroResultsInRoi(Roi roi, Objects3DIntPopulation microPop, ImagePlus imgMicro, ImageHandler imhMicro, 
                                       AtomicInteger microLabel, ImagePlus imgVesselMaskRoiDil, ImagePlus imgVesselSkel, 
                                       ImageFloat vesselDistMap, ImageFloat vesselDistMapInv, Calibration cal, double vesselVol, 
                                       BufferedWriter microResults, BufferedWriter globalResults, 
                                       String imgName, String roiName) throws IOException {
        
        System.out.println("Computing microglia parameters...");
        
        // Get microglia in ROI
        Objects3DIntPopulation microPopInRoi = getPopInsideRoi(microPop, imgMicro, roi);
        
        // Get skeleton in (translated) dilated ROI
        ImagePlus imgVesselSkelRoiDil = clearOutsideRoi(imgVesselSkel, roi, true, cal);
        
        int nbVAM=0, nbVTM=0, nbVDM = 0;
        for (Object3DInt micro: microPopInRoi.getObjects3DInt()) {
            double microVol = new MeasureVolume(micro).getVolumeUnit();
            microResults.write(imgName+"\t"+roiName+"\t"+microLabel.get()+"\t"+microVol);
            micro.setLabel(microLabel.getAndIncrement());
            if(vesselVol == 0) {
                microResults.write("\t"+Double.NaN+"\t"+Double.NaN+"\t"+Double.NaN+"\t"+Double.NaN+"\n");
                microResults.flush();
                nbVDM++;
            } else {  
                Object3DInt vesselObj = new Object3DInt(ImageHandler.wrap(imgVesselMaskRoiDil));   
                double colocVol = getColocVol(micro, vesselObj, cal);
                
                double centroidDist = vesselDistMapInv.getPixel(new MeasureCentroid​(micro).getCentroidAsPoint());
                double borderDist = vesselDistMapInv.getPixel(new Measure2Distance(micro, vesselObj).getBorder1Pix());
                
                Object3DInt skelObj = new Object3DInt(ImageHandler.wrap(imgVesselSkelRoiDil));
                double vesselDiam = 2*vesselDistMap.getPixel(new Measure2Distance(micro, skelObj).getBorder2Pix());
                
                microResults.write("\t"+colocVol+"\t"+centroidDist+"\t"+borderDist+"\t"+vesselDiam+"\n");
                microResults.flush();

                if(colocVol == 0) nbVDM++;
                else if(colocVol != 0 && centroidDist == 0) nbVAM++;
                else if(colocVol != 0 && centroidDist != 0) nbVTM++;
            }            
        }
        
        microPopInRoi.drawInImage(imhMicro);
                
        globalResults.write("\t"+microPopInRoi.getNbObjects()+"\t"+nbVAM+"\t"+nbVTM+"\t"+nbVDM);
        globalResults.flush();
        
        closeImage(imgVesselSkelRoiDil);
    }
    
    
    private void saveEndoResultsInRoi(Roi roi, Objects3DIntPopulation endoPop, ImagePlus imgEndo, ImageHandler imhEndo, double roiVol, 
            double vesselTotalLength, BufferedWriter globalResults) throws IOException {
        
        System.out.println("Computing endothelial nuclei parameters...");
        
        // Get endothelial nuclei in translated ROI
        Objects3DIntPopulation endoPopInRoi = getPopInsideRoi(endoPop, imgEndo, roi);
        endoPopInRoi.drawInImage(imhEndo);

        // Write endothelial nuclei parameters in global results file
        int nbEndoInRoi = endoPopInRoi.getNbObjects();
        if(vesselTotalLength != 0)
            globalResults.write("\t"+nbEndoInRoi+"\t"+nbEndoInRoi/roiVol*1e6+"\t"+nbEndoInRoi/vesselTotalLength);
        else
            globalResults.write("\t"+nbEndoInRoi+"\t"+nbEndoInRoi/roiVol*1e6+"\t"+Double.NaN);
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
    
    
    private void saveImages(ImageStack stackTagSkel, ImageHandler imhVessels, ImageHandler imhMicro, ImageHandler imhEndo, 
            ImagePlus imgVessels, ImagePlus imgMicro, ImagePlus imgEndo, Calibration cal, String rootName, String outDir)  {
        
        ImagePlus imgTagSkel = new ImagePlus("", stackTagSkel);
        IJ.run(imgTagSkel, "Fire", null);
        imgTagSkel.resetDisplayRange();
        ImagePlus[] imgStack1 = {imgTagSkel, null, null, imhVessels.getImagePlus()};
        ImagePlus imgMerge1 = new RGBStackMerge().mergeHyperstacks(imgStack1, true);
        imgMerge1.setCalibration(cal);
        new FileSaver(imgMerge1).saveAsTiff(outDir+rootName+"_skeleton.tif");
        closeImage(imgTagSkel);
        closeImage(imgMerge1);
        
        ImagePlus[] imgStack2 = {imhVessels.getImagePlus(), null, null, imgVessels, null, null};
        if (imgMicro != null) {
            imgStack2[1] = imhMicro.getImagePlus();
            imgStack2[4] = imgMicro;
        }
        if (imgEndo != null) {
            imgStack2[2] = imhEndo.getImagePlus();
            imgStack2[5] = imgEndo;
        }
        ImagePlus imgMerge2 = new RGBStackMerge().mergeHyperstacks(imgStack2, true);
        imgMerge2.setC(1);
        imgMerge2.setDisplayRange(0, 1);
        imgMerge2.setC(2);
        imgMerge2.setDisplayRange(0, 1);
        imgMerge2.setC(4);
        IJ.run(imgMerge2, "Grays", "");
        imgMerge2.setC(5);
        IJ.run(imgMerge2, "Grays", "");
        imgMerge2.setCalibration(cal);
        new FileSaver(imgMerge2).saveAsTiff(outDir+rootName+".tif");
        closeImage(imgMerge2);
    }
    
    
}
