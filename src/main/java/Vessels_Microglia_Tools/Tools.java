package Vessels_Microglia_Tools;


import Vessels_Microglia_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Vessels_Microglia_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.util.ThreadUtil;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.stream.IntStream;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
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
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.SkeletonResult;
import trainableSegmentation.WekaSegmentation;
import weka.core.Utils;



public class Tools {
        
    // min size for objects
    public double minVessels = 100;
    public double minMicro = 50;
    // max size for objects
    public double maxVessels = Double.MAX_VALUE;
    public double maxMicro = Double.MAX_VALUE;
    public Calibration cal = new Calibration();
    public double pixVol;
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    public String modelMicro = "cyto2";
    public int cellposeDiamMicro = 30;
    public double cellposeStitchTh = 0.5;
    public String modelVessel = "vessels";
     public int cellposeDiamVessel = 40;
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
     /**
     * Find images extension
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
               case "nd" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }

        
   /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + f);
        }
        Collections.sort(images);
        return(images);
    }
    
      /**
     * Find weka model in folder
     * @param imagesFolder
     * @param model
     * @return 
     */
    public String findWekaModel(String imagesFolder, String model) {
        String modelFile = "";
        ArrayList<String> files = findImages(imagesFolder, "model");
        if (files.isEmpty()) {
            System.out.println("No Model found in "+imagesFolder);
            return null;
        }
        for (String f : files) {
            String name = FilenameUtils.getBaseName(f);
            if (model.equals("vessel") && name.equals("vessel")) {
                modelFile = f;
                break;
            }
            if (model.equals("microglia") && name.equals("microglia")) {
                modelFile = f;
                break;
            }
        }        
        return modelFile;
    }
        
    
    /**
     * Dialog
     */
    public String[] dialog(String imagesDir, String[] channels) {
        String[] chNames = {"Vessels : ", "Microglia : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 15, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[n]);
        }
        gd.addDirectoryField​("Cellpose env directory", cellposeEnvDir);
        gd.addMessage("Models for vessels and microglia detections", Font.getFont("Monospace"), Color.blue);
        gd.addFileField("Model for vessels", modelVessel);
        gd.addFileField("Model for microglia", modelMicro);
        gd.addMessage("Size filter", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min vessels size (µm3) : ", minVessels, 2);
        gd.addNumericField("Min microglia size (µm3) : ", minMicro, 2);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("xy calibration (µm)", cal.pixelHeight,3);
        gd.addNumericField("z calibration (µm)", cal.pixelDepth,3);
        
        gd.showDialog();
        String[] chChoices = new String[chNames.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();
        
        if (gd.wasCanceled())
                chChoices = null;
        if (modelVessel.isEmpty() || modelMicro.isEmpty()) {
            modelVessel = gd.getNextString();
            modelMicro = gd.getNextString();
        }
        minVessels = gd.getNextNumber();
        minMicro = gd.getNextNumber();
        cal.pixelHeight = gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight;
        cal.pixelDepth = gd.getNextNumber();
        pixVol = cal.pixelHeight*cal.pixelWidth*cal.pixelDepth;
        return(chChoices);
    } 
    
     /**
     * Find channels name and None to end of list
     * @param imageName
     * @param meta
     * @param reader
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
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
    

    public void setCalibration(ImagePlus imp) {
        imp.setCalibration(cal);
    }
     
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
    
      public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
      
      /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.setVoxelSizeXY(cal.pixelWidth);
        pop.setVoxelSizeZ(cal.pixelDepth);
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
      
    /**
     * 3D Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus median3D_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
        clij2.release(imgCLMed);
       return(imgMed);
    }
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Save images
     * @param vesselsObj
     * @param microPop
     * @param imgs
     * @return 
     */
    public ImagePlus[] createImageObjects(Object3DInt vesselsObj, Objects3DIntPopulation microPop, ImagePlus[] imgs, ImagePlus imgVessels) {
        ImageHandler imhVessels = (imgs[0] == null) ? ImageHandler.wrap(imgVessels).createSameDimensions() : ImageHandler.wrap(imgs[0]);
        ImageHandler imhMicro = (imgs[1] == null) ? imhVessels.createSameDimensions() : ImageHandler.wrap(imgs[1]);
        vesselsObj.drawObject(imhVessels, 255);
        microPop.drawInImage(imhMicro);
        ImagePlus[] imgObjects = {imhVessels.getImagePlus(), imhMicro.getImagePlus()};
        imhVessels.closeImagePlus();
        imhMicro.closeImagePlus();
        return(imgObjects);
    }
    
    
     /**
     * compute local thickness
     * @param img
     * @param inverse
     * @return imgMap
    **/
    public ImageFloat localThickness3D (ImagePlus img, boolean inverse) {
        IJ.showStatus("Computing distance map...");
        img.setCalibration(cal);
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0, inverse, ThreadUtil.getNbCpus());
        return(edt);
    }
    
    /**
     * Clij2 skeletonize 3D
     * @param img
     * @return 
     */
    public ImagePlus skeletonize3D(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLSkel = clij2.create(imgCL);
        BoneJSkeletonize3D skel = new BoneJSkeletonize3D();
        skel.bonejSkeletonize3D(clij2,imgCL, imgCLSkel);
        clij2.release(imgCL);
        ImagePlus imgSkel = clij2.pull(imgCLSkel);
        clij2.release(imgCLSkel);
        return(imgSkel);
    }
    
    /**
     * Reset labels of the objects composing a population
     * @param pop
     */
    public void resetLabels(Objects3DPopulation pop) {
        for(int i=0; i < pop.getNbObjects(); i++) {
            pop.getObject(i).setValue(i+1);
        }
    }
    
    public void preprocessFile(String imageDir, String processDir, ArrayList<String> imageFiles, String[] chsName, String[] channels) throws Exception {
        try {
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                ImageProcessorReader reader = new ImageProcessorReader();
                reader.setId(f);
                ImporterOptions options = new ImporterOptions();
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setId(f);
                options.setSplitChannels(true);
                int indexCh = ArrayUtils.indexOf(chsName, channels[0]);
                // Open vessel channel
                options.setCBegin(0, indexCh);
                options.setCEnd(0, indexCh);
                ImagePlus imgVessels = BF.openImagePlus(options)[0];
                setCalibration(imgVessels);
                IJ.saveAs(imgVessels, "Tiff", processDir+rootName+"-Vessels.tif");
                closeImages(imgVessels);

                // Open microglia channel
                indexCh = ArrayUtils.indexOf(chsName, channels[1]);
                options.setCBegin(0, indexCh);
                options.setCEnd(0, indexCh);
                ImagePlus imgMicro = BF.openImagePlus(options)[0];
                setCalibration(imgMicro);
                IJ.saveAs(imgMicro, "Tiff", processDir+rootName+"-Micro.tif");
                closeImages(imgMicro);
            }
        }
        catch (Exception e) { throw e; }
    }
    
     /*
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus imgIn, String model, int diameter) throws IOException{
       ImagePlus img = new Duplicator().run(imgIn);
       
       // Define CellPose settings
       CellposeTaskSettings settings = new CellposeTaskSettings(model, 1, diameter, cellposeEnvDir);
       settings.setStitchThreshold(cellposeStitchTh);
       settings.useGpu(true);
       
       // Run CellPose
       CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
       ImagePlus imgOut = cellpose.run();
       imgOut = imgOut.resize(imgIn.getWidth(), imgIn.getHeight(), "none");
       imgOut.setCalibration(cal);
       closeImages(img);
       // remove small objects
        if (model.equals(modelVessel))
            IJ.run(imgOut, "Options...", "iterations=1 count=1 do=Erode stack");
      
        // Get objects population
        Objects3DIntPopulation pop = getPopFromImage(imgOut);
        System.out.println(pop.getNbObjects()+" before size filtering");
        if (model.equals(modelVessel))
            popFilterSize(pop, minVessels, maxVessels);
        else
            popFilterSize(pop, minMicro, maxMicro);
        closeImages(imgOut);
        System.out.println(pop.getNbObjects()+" after size filtering");
        return(pop);
    }
     
    /**
     * Find vessels length and mean diameter
     */ 
    private HashMap skeletonParameters(ImagePlus vesselsSkelImg, ImagePlus vesselimgMap) {
        System.out.println("Computing vessel parameters ...");
        HashMap<String, Double> params = new HashMap<>();
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",vesselsSkelImg);
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
        double[] branchLengths = skeletonResults.getAverageBranchLength();
        int[] branchNumbers = skeletonResults.getBranches();
        int totalLength = 0;
        for (int i = 0; i < branchNumbers.length; i++)
            totalLength += branchNumbers[i] * branchLengths[i];
        params.put("totalLength", (double)totalLength);
        params.put("meanLength", StatUtils.mean(branchLengths));
        params.put("lengthLongestBranch", StatUtils.max(skeletonResults.getMaximumBranchLength()));
        params.put("nbBranches", (double)IntStream.of(skeletonResults.getBranches()).sum());
        params.put("nbJunctions", (double)IntStream.of(skeletonResults.getJunctions()).sum());
        params.put("nbEndpoints", (double)IntStream.of(skeletonResults.getEndPoints()).sum());
        
        ArrayList<Point> AllVoxels = skeletonResults.getListOfSlabVoxels();
        ArrayList<Point> JunctionVoxels = skeletonResults.getListOfJunctionVoxels();
        AllVoxels.addAll(JunctionVoxels);
        double MeanDiam = 0;
        for (Point pt : AllVoxels) {
            vesselimgMap.setZ(pt.z);
            MeanDiam += vesselimgMap.getProcessor().getPixelValue(pt.x, pt.y);
        }
        params.put("meanDiameter", MeanDiam/AllVoxels.size());
        return(params);
    } 
    
    /**
     * remove object outside roi
     */
     private Object3DInt removeObjOutsideRoi(Objects3DIntPopulation pop, ImagePlus img, Roi roi) {
        // get pop as one 3D Object
        ImageHandler imhVessels = ImageHandler.wrap(img).createSameDimensions();
        for (Object3DInt obj : pop.getObjects3DInt())
            obj.drawObject(imhVessels, 255);
        // take only inside roi
        ImagePlus imgMask = imhVessels.getImagePlus();
        imgMask.setRoi(roi);
        IJ.run(imgMask, "Clear Outside", "stack");
        imgMask.deleteRoi();
        imgMask.setCalibration(cal);
        Object3DInt vesselObj = new Object3DInt(ImageHandler.wrap(imgMask)); 
        return(vesselObj);
     }
     
     /**
     * remove objects outside roi
     */
     private Objects3DIntPopulation removePopOutsideRoi(Objects3DIntPopulation pop, ImagePlus img, Roi roi) {
        ImageHandler imhVessels = ImageHandler.wrap(img).createSameDimensions();
        pop.drawInImage(imhVessels);
        // take only inside roi
        ImagePlus imgMask = imhVessels.getImagePlus();
        imgMask.setRoi(roi);
        IJ.run(imgMask, "Clear Outside", "stack");
        imgMask.deleteRoi();
        imgMask.setCalibration(cal);
        Objects3DIntPopulation microPop = getPopFromImage(imgMask); 
        microPop.resetLabels();
        return(microPop);
     }
     
     /*
     Get coloc volume of microglial cell and vessel 
     */
     private double getColocVol(Object3DInt microObj, Object3DInt vessel) {
        Measure2Colocalisation coloc = new Measure2Colocalisation(microObj, vessel);
        double colVol = coloc.getValue(Measure2Colocalisation.COLOC_VOLUME);
        return(colVol*pixVol);
     }
   
    /** compute parameters and save images objects 
     * 
     * @param vesselsPop
     * @param microPop
     * @param roi
     * @param imgVessels
     * @param imgMicro
     * @param rootName
     * @param roiName
     * @param outPutResults
     * @return 
     * @throws java.io.IOException 
     */
    public ImagePlus[] compute_param(Objects3DIntPopulation vesselsPop, Objects3DIntPopulation microPop, Roi roi,
            ImagePlus imgVessels, ImagePlus imgMicro, ImagePlus[] imgs, String rootName, String outDir,
            BufferedWriter outPutResults) throws IOException {
        // get vessels pop as one 3D Object in roi
        Object3DInt vesselInRoiObj = removeObjOutsideRoi(vesselsPop, imgVessels, roi);
        double vesselsVol = new MeasureVolume(vesselInRoiObj).getVolumeUnit();
        // Compute distance map
        ImageHandler imhMap = ImageHandler.wrap(imgVessels).createSameDimensions();
        vesselInRoiObj.drawObject(imhMap, 255);
        ImagePlus imgMask = imhMap.getImagePlus();
        imhMap.closeImagePlus();
        ImageFloat vesselimgMap = localThickness3D(imgMask, false);
        ImageFloat vesselimgMapInv = localThickness3D(imgMask, true);
        System.out.println("compute vessels skeleton "+roi.getName()+" ...");
        ImagePlus vesselsSkelImg = skeletonize3D(imgMask);
        Object3DInt vesselsSkelObj = new Object3DInt(ImageHandler.wrap(vesselsSkelImg));
        closeImages(imgMask);
        IJ.run(vesselsSkelImg, "8-bit","");
        HashMap skelParams = skeletonParameters(vesselsSkelImg, vesselimgMap.getImagePlus());
        // get microglia pop in roi
        Objects3DIntPopulation microRoiPop = removePopOutsideRoi(microPop, imgMicro, roi);  
        int index = 0;
        for (Object3DInt microObj : microRoiPop.getObjects3DInt()) {
            double microVol = new MeasureVolume(microObj).getVolumeUnit();
            Point3D pt = new MeasureCentroid​(microObj).getCentroidAsPoint();
            double CenterDist = vesselimgMapInv.getPixel(pt);
            double BorderDist = new Measure2Distance(microObj, vesselInRoiObj).getValue(Measure2Distance.DIST_BB_UNIT);
            VoxelInt voxelBorder = new Measure2Distance(microObj, vesselsSkelObj).getBorder2Pix();
            double radius = vesselimgMap.getPixel(voxelBorder);
            double microColocVol = getColocVol(microObj, vesselInRoiObj);
            if (index == 0)
                outPutResults.write(rootName+"\t"+roi.getName()+"\t"+microObj.getLabel()+"\t"+microVol+"\t"+microColocVol+"\t"+CenterDist+"\t"+BorderDist+"\t"+radius+"\t"+vesselsVol+"\t"+
                    skelParams.get("totalLength")+"\t"+skelParams.get("meanLength")+"\t"+skelParams.get("lengthLongestBranch")+"\t"+
                    skelParams.get("nbBranches")+"\t"+skelParams.get("nbJunctions")+"\t"+skelParams.get("nbEndpoints")+"\t"+skelParams.get("meanDiameter")+"\n");
            else
                outPutResults.write("\t\t"+microObj.getLabel()+"\t"+microVol+"\t"+microColocVol+"\t"+CenterDist+"\t"+BorderDist+"\t"+radius+"\t\t\t\t\t\t\t\t\n");
            outPutResults.flush();  
            index++;
        }
        vesselimgMap.closeImagePlus();
        vesselimgMapInv.closeImagePlus();
        // Draw image objects
        System.out.println("Create image objects");
        imgs = createImageObjects(vesselInRoiObj, microRoiPop, imgs, imgVessels);
        return(imgs);
    }      
}
