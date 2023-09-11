import Vessels_Cellpose_Microglia_Tools.QuantileBasedNormalization;
import Vessels_Cellpose_Microglia_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
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

/**
 * Detect vessels and microglial cells with Cellpose
 * Analyze vessels skeleton
 * Compute distance between each microglial cell and its nearest vessel
 * @authors Philippe Mailly & Héloïse Monnet
 */
public class Vessels_Cellpose_Microglia implements PlugIn {

    private Vessels_Cellpose_Microglia_Tools.Tools tools = new Tools();
    
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            }
            
            String imageDir = IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                return;
            }
            
            // Find images with fileExt extension
            String fileExt = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] channelNames = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] channels = tools.dialog(imageDir, channelNames);
            if(channels == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            // Create output folder for results files and images
            String outDir = imageDir + File.separator + "Results_" + tools.cellposeModelVessel + "_" + 
                    new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write headers in results files
            FileWriter fwMicroResults = new FileWriter(outDir + "microgliaResults.csv", false);
            BufferedWriter microResults = new BufferedWriter(fwMicroResults);
            microResults.write("Image name\tRoi name\tCell ID\tCell volume (µm3)\tCell volume coloc with vessel (µm3)"
                    + "\tCell centroid distance to closest vessel (µm)\tCell border distance to closest vessel (µm)"
                    + "\tClosest vessel diameter (µm)\n");
            microResults.flush();
            FileWriter fwGlobalResults = new FileWriter(outDir + "globalResults.csv", false);
            BufferedWriter globalResults = new BufferedWriter(fwGlobalResults);
            globalResults.write("Image name\tRoi name\tRoi volume (µm3)\tVessels total volume (µm3)\tVessels density\tVessels total length (µm)"
                    + "\tVessels mean length (µm)\tLongest branch length (µm)\tNb branches\tNb junctions\tNb endpoints"
                    + "\tVessels mean diameter (µm)\tVessels std diameter (µm)\tVessels min diameter (µm)\tVessels max diameter (µm)"
                    + "\tNb cells\tNb VAM\tNb VTM\tNb VDM\n");
            globalResults.flush();
            
            // Create output folder for preprocessed images
            String processDir = imageDir + File.separator + "Preprocessed" + File.separator;
            if (!Files.exists(Paths.get(processDir))) {
                new File(processDir).mkdir();
                
                // Preprocess images
                tools.print("--- NORMALIZING IMAGES ---");
                tools.preprocessFiles(imageFiles, processDir, channelNames, channels);
                // Normalize vessels
                QuantileBasedNormalization qbn = new QuantileBasedNormalization();
                qbn.run(processDir, imageFiles, "-Vessels");
                // Normalize microglia (if channel exists)
                if (!channels[2].equals("None") && tools.microMethod.equals("Cellpose"))
                    qbn.run(processDir, imageFiles, "-Micro");
                tools.print("Normalisation done");
            }
            
            IJ.setForegroundColor(255, 255, 255);
            IJ.setBackgroundColor(0, 0, 0);
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ---");
                
                // Read vessels and microglia channels
                tools.print("- Reading channels -");
                ImagePlus imgVessels = IJ.openImage(processDir + rootName + "-Vessels-normalized.tif");
                ImagePlus imgMicro = null;
                if (!channels[2].equals("None")) {
                    String fileName = (tools.microMethod.equals("Cellpose")) ? processDir+rootName+"-Micro-normalized.tif" : processDir+rootName+"-Micro.tif"; 
                    imgMicro = IJ.openImage(fileName);
                }
                
                // Load ROIs (if provided)
                tools.print("- Loading ROIs -");
                List<Roi> rois = new ArrayList();
                String roiName = imageDir + File.separator + rootName; 
                roiName = new File(roiName + ".zip").exists() ? roiName + ".zip" : roiName + ".roi";
                if (new File(roiName).exists()) {
                    RoiManager rm = new RoiManager(false);
                    rm.runCommand("Open", roiName);
                    rois = Arrays.asList(rm.getRoisAsArray());
                } else {
                    Roi roi = new Roi(0, 0, imgVessels.getWidth() - 1 , imgVessels.getHeight() - 1);
                    roi.setName("whole image");
                    rois.add(roi);
                }
                
                // Segment vessels with Cellpose
                tools.print("- Segmenting vessels -");
                Objects3DIntPopulation vesselPop = tools.cellposeDetection(imgVessels, tools.cellposeModelsPath+tools.cellposeModelVessel, tools.cellposeDiamVessel, tools.minVesselVol, tools.maxVesselVol);
                
                // Segment microglia with Cellpose or Otsu thresholding
                Objects3DIntPopulation microPop = new Objects3DIntPopulation();
                if (imgMicro != null) {
                    tools.print("- Segmenting microglial cells -");
                    microPop = (tools.microMethod.equals("Cellpose")) ? tools.cellposeDetection(imgMicro, tools.cellposeModelMicro, tools.cellposeDiamMicro, tools.minMicroVol, tools.maxMicroVol) :
                        tools.microgliaSegmentation(imgMicro, tools.minMicroVol, tools.maxMicroVol);
                }
               
                // Save results
                tools.saveResults(vesselPop, microPop, rois, imgVessels, imgMicro, microResults, globalResults, rootName, outDir);
                
                tools.closeImage(imgVessels);
                if (imgMicro != null) 
                    tools.closeImage(imgMicro);
            }
            microResults.close();
            globalResults.close();
        } catch (IOException ex) {
            Logger.getLogger(Vessels_Cellpose_Microglia.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(Vessels_Cellpose_Microglia.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("All done!");
    }
}
