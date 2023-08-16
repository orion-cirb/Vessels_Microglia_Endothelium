/*
 * Analyze microglia on vessels
 * Use cellpose with vessel model
 * Author Philippe Mailly / Heloïse Monnet
 */


import Vessels_Microglia_Tools.QuantileBasedNormalization;
import Vessels_Microglia_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import loci.common.services.ServiceFactory;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;


public class Vessels_Cellpose_Microglia implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    private Vessels_Microglia_Tools.Tools tools = new Tools();
    
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            String imageDir = IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            // Find images with fileExt extension
            String fileExt = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find chanels, image calibration
            reader.setId(imageFiles.get(0));
            String[] channelNames = tools.findChannels(imageFiles.get(0), meta, reader);
            tools.findImageCalib(meta);
            String[] channels = tools.dialog(imageDir, channelNames);
            if(channels == null)
                return;
            
            // Create output folder
            DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
            Date date = new Date();
            String outDirResults = imageDir + File.separator+ "Results_" + tools.modelVessel+"_"+dateFormat.format(date) +  File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter results = new FileWriter(outDirResults + "Results.xls", false);
            BufferedWriter outPutResults = new BufferedWriter(results);
            outPutResults.write("ImageName\tRoi name\tRoi volume (µm3)\t#Microglia cell\tCell volume (µm3)\tCell volume coloc with vessel(µm3)\tCentroid cell distance to vessel (µm)\tBorder cell distance to vessel (µm)\tVessel diameter (µm)\tVessel volume (µm3)"
                    + "\tTotal vessel length (µm)\tMean length (µm)\tLength longest Branch (µm)\tNb branches\tNb junctions\tnb endpoints\tMean vessel diameter (µm)"
                    + "\tStd vessel diameter (µm)\tMin vessel diameter (µm)\tMax vessel diameter (µm)\n");
            outPutResults.flush();
            
           // Create output folder for preprocessed file
            String processDir = inDir + File.separator + "Preprocessed"+ File.separator;
            File procDir = new File(processDir);
            if (!Files.exists(Paths.get(processDir))) {
                procDir.mkdir();
                // Normalise all images
                IJ.showStatus("Normalisation starting...");
                tools.preprocessFile(imageDir, processDir, imageFiles, channelNames, channels);
                QuantileBasedNormalization qbn = new QuantileBasedNormalization();
                // vessels normalization
                qbn.run(processDir, imageFiles, "-Vessels");
                // microglia if channel exists
                if (tools.cellsDetection.equals("CellPose") && !channels[2].equals("None"))
                    qbn.run(processDir, imageFiles, "-Micro");
                IJ.showStatus("Normalisation done");
            }
               
                    
            // Do cellpose detection on all files
            IJ.setForegroundColor(255, 255, 255);
            IJ.setBackgroundColor(0, 0, 0);

            System.out.println("Segmentations ...");
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                String fileName = processDir+rootName+"-Vessels-normalized.tif"; 
                String roiFileName = null;
                if (new File(imageDir+rootName+".zip").exists())
                    roiFileName = imageDir+rootName+".zip";
                else if (new File(imageDir+rootName+".roi").exists())
                    roiFileName = imageDir+rootName+".roi";
                reader.setId(fileName);
                ImagePlus imgVessels = IJ.openImage(fileName);
                ImagePlus imgMicro = null;
                if (!channels[2].equals("None")) {
                    fileName = (tools.cellsDetection.equals("CellPose")) ? processDir+rootName+"-Micro-normalized.tif" : processDir+rootName+"-Micro.tif"; 
                    imgMicro = IJ.openImage(fileName);
                }
                
                //do vessel segmentation
                System.out.println("Vessels segmentation starting image " + rootName +" ...");
                Objects3DIntPopulation vesselPop = tools.cellposeDetection(imgVessels, tools.modelVessel, 35);
                
                //do microglia segmentation
                Objects3DIntPopulation microPop = new Objects3DIntPopulation();
                if (imgMicro != null) {
                    System.out.println("Microglia segmentation starting image " + rootName + " ...");
                    microPop = (tools.cellsDetection.equals("CellPose")) ? tools.cellposeDetection(imgMicro, tools.modelMicro, 30) :
                        tools.microgliaCells(imgMicro);
                }
               
                // Load Roi
                List<Roi> rois = new ArrayList();
                if (roiFileName == null) {
                    Roi roi = new Roi(0, 0, imgVessels.getWidth() - 1 , imgVessels.getHeight() - 1);
                    roi.setName("whole image");
                    rois.add(roi);
                }
                else {
                    RoiManager rm = new RoiManager(false);
                    rm.runCommand("Open", roiFileName);
                    rois = Arrays.asList(rm.getRoisAsArray());
                }
                ImagePlus[] imgs = {null, null};
               
                for (Roi roi : rois) {
                    // Compute parameters
                    System.out.println("computing parameters ...");
                    imgs = tools.compute_param(vesselPop, microPop, roi, imgVessels, imgMicro, imgs, rootName, outDirResults, outPutResults);
                } 
                IJ.run(imgs[0], "Red", "");
                IJ.run(imgs[1], "Green", "");
                ImagePlus[] imgColors = {imgs[0], imgs[1], null, imgVessels, imgMicro};
                ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
                imgObjects.setCalibration(tools.cal);
                FileSaver imgObjectsFile = new FileSaver(imgObjects);
                imgObjectsFile.saveAsTiff(outDirResults+rootName+"_cells.tif");
                tools.closeImages(imgObjects);
                tools.closeImages(imgVessels);
                if (imgMicro != null) 
                        tools.closeImages(imgMicro);
            }

            outPutResults.close();

        } catch (IOException ex) {
            Logger.getLogger(Vessels_Cellpose_Microglia.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(Vessels_Cellpose_Microglia.class.getName()).log(Level.SEVERE, null, ex);
        }
            
            IJ.showStatus("Process done");
    }
}
