import Vessels_Cellpose_Microglia_Tools.QuantileBasedNormalization;
import Vessels_Cellpose_Microglia_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
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
import java.util.Date;
import java.util.List;
import loci.common.services.ServiceFactory;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

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
            
            // Find images with tif extension
            ArrayList<String> imageFiles = tools.findImages(imageDir, "tif");
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with tif extension");
                return;
            }
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            
            // Find channel names
            String[] channelNames = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] channels = tools.dialog(imageDir, channelNames);
            if(channels == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            // Create output folder for results files and images
            String vesselMethodName = (tools.vesselSegMethod == "Cellpose")? tools.cellposeModelVessel : tools.vesselThMethod;
            String microMethodName = (channels[1].equals("None")) ? "" : tools.microThMethod + "_";
            String outDir = imageDir + File.separator + "Results_" + vesselMethodName + "_" + microMethodName +
                    new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write headers in results files
            FileWriter fwMicroResults = new FileWriter(outDir + "microgliaResults.csv", false);
            BufferedWriter microResults = new BufferedWriter(fwMicroResults);
            microResults.write("Image name\tROI name\tCell ID\tCell volume (µm3)\tCell volume coloc with vessel (µm3)"
                    + "\tCell centroid distance to closest vessel (µm)\tCell border distance to closest vessel (µm)"
                    + "\tClosest vessel diameter (µm)\n");
            microResults.flush();
            FileWriter fwGlobalResults = new FileWriter(outDir + "globalResults.csv", false);
            BufferedWriter globalResults = new BufferedWriter(fwGlobalResults);
            globalResults.write("Image name\tImage XY calibration (µm)\tROI name\tROI translation\tROI volume (µm3)\tVessels total volume (µm3)\tVessels density\tVessels total length (µm)"
                    + "\tVessels mean length (µm)\tLongest branch length (µm)\tNb branches\tNb junctions\tNb endpoints"
                    + "\tVessels mean diameter (µm)\tVessels std diameter (µm)\tVessels min diameter (µm)\tVessels max diameter (µm)"
                    + "\tNb cells\tNb VAM\tNb VTM\tNb VDM\n");
            globalResults.flush();
            
            // Preprocess vessels channel
            String processDir = imageDir + File.separator + "Preprocessed" + File.separator;
            if (!Files.exists(Paths.get(processDir))) {
                tools.print("--- NORMALIZING IMAGES ---");
                // Create output folder for preprocessed files
                new File(processDir).mkdir();
                // Save vessels channel
                tools.preprocessFiles(reader, imageFiles, processDir, channelNames, channels);
                // Normalize vessels channel
                QuantileBasedNormalization qbn = new QuantileBasedNormalization();
                qbn.run(processDir, imageFiles, "-Vessels");
                tools.print("Normalisation done");
            }
            
            IJ.setForegroundColor(255, 255, 255);
            IJ.setBackgroundColor(0, 0, 0);
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ---");
                reader.setId(f);
                
                // Find image calibration
                Calibration cal = tools.findImageCalib(meta);
                                
                // Open vessels and microglia channels
                tools.print("- Opening channels -");
                ImagePlus imgVessels = IJ.openImage(processDir + rootName + "-Vessels-normalized.tif");
                ImagePlus imgMicro = null;
                if (!channels[1].equals("None")) {
                    ImporterOptions options = new ImporterOptions();
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setQuiet(true);
                    
                    int indexCh = ArrayUtils.indexOf(channelNames, channels[1]);
                    imgMicro = BF.openImagePlus(options)[indexCh];
                }
                
                // Load ROIs (if provided)
                tools.print("- Loading ROIs -");
                List<Roi> rois = tools.loadRois(imageDir + File.separator + rootName, imgVessels, rootName);
                
                // Segment vessels
                tools.print("- Segmenting vessels -");
                ImagePlus imgVesselMask = tools.vesselSegmentation(imgVessels, cal);
                
                // Segment microglia
                Objects3DIntPopulation microPop = new Objects3DIntPopulation();
                if (imgMicro != null) {
                    tools.print("- Segmenting microglia -");
                    microPop = tools.cellsSegmentation(imgMicro, cal);
                }
               
                // Save results
                tools.saveResults(imgVesselMask, microPop, rois, imgVessels, imgMicro, cal, microResults, globalResults, rootName, outDir);

                tools.closeImage(imgVessels);
                if (imgMicro != null) tools.closeImage(imgMicro);
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
