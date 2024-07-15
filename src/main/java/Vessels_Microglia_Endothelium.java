import Vessels_Microglia_Endothelium_Tools.QuantileBasedNormalization;
import Vessels_Microglia_Endothelium_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.LUT;
import java.awt.Color;
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
import mcib3d.image3d.ImageFloat;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/**
 * Detect vessels, microglial cells and endothelial nuclei
 * Analyze vessels structure
 * Compute distance between microglial cells and their nearest vessel
 * Give number of endothelial nuclei
 * @authors Héloïse Monnet & Philippe Mailly
 */
public class Vessels_Microglia_Endothelium implements PlugIn {

    private Vessels_Microglia_Endothelium_Tools.Tools tools = new Tools();
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            }
            
            String imageDir = IJ.getDirectory("Choose images directory");
            if (imageDir == null) {
                return;
            }
            
            // Find images with tif extension
            ArrayList<String> imageFiles = tools.findImages(imageDir, "tif");
            if (imageFiles.isEmpty()) {
                IJ.showMessage("ERROR", "No image found with tif extension in " + imageDir + " folder");
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
                return;
            } else if(channels[0] == "None") {
                IJ.showMessage("ERROR", "Vessels channel not defined");
                return;
            }
            
            // Create output folder for results files and images
            String vesselNorm = tools.vesselNormalization? "Norm_" : "";
            String vesselMethodName = (tools.vesselSegMethod == "Cellpose")? tools.cellposeModelVessel : tools.vesselThMethod;
            String microMethodName = (channels[1].equals("None")) ? "" : tools.microThMethod + "_";
            String outDir = imageDir + File.separator + "Results_" + vesselNorm + vesselMethodName + "_" + microMethodName +
                            new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write headers in results files
            // Global
            FileWriter fwGlobalResults = new FileWriter(outDir + "globalResults.csv", false);
            BufferedWriter globalResults = new BufferedWriter(fwGlobalResults);
            // Vessels
            FileWriter fwVesselResults = new FileWriter(outDir + "vesselsResults.csv", false);
            BufferedWriter vesselResults = new BufferedWriter(fwVesselResults);
            // Microglia
            BufferedWriter microResults = null;
            if(channels[1] != "None") {
                FileWriter fwMicroResults = new FileWriter(outDir + "microgliaResults.csv", false);
                microResults = new BufferedWriter(fwMicroResults);
            }
            tools.writeHeaders(channels, globalResults, vesselResults, microResults);
            
            // If asked in dialog box, normalize vessels channel
            String normDir = imageDir + File.separator + "Normalization" + File.separator;
            if (tools.vesselNormalization && !Files.exists(Paths.get(normDir))) {
                tools.print("--- NORMALIZING IMAGES ---");
                // Create output folder for normalized files
                new File(normDir).mkdir();
                // Save images vessels channel
                tools.saveChannel(normDir, imageFiles, ArrayUtils.indexOf(channelNames, channels[0]), "-vessels.tif");
                // Normalize images vessels channel
                new QuantileBasedNormalization().run(normDir, imageFiles, "-vessels");
                // Delete images vessels channel
                tools.deleteChannel(normDir, imageFiles, "-vessels.tif");
                tools.print("Normalization done");
            }
            
            IJ.setForegroundColor(255, 255, 255);
            IJ.setBackgroundColor(0, 0, 0);
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ---");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setQuiet(true);
                
                // Find image calibration
                Calibration cal = tools.findImageCalib(meta);
                                
                // Open channels
                tools.print("- Opening channels -");
                LUT lut = LUT.createLutFromColor(Color.gray);
                
                ImagePlus imgVessels = null;
                if (tools.vesselNormalization) { 
                    imgVessels = IJ.openImage(normDir + rootName + "-vessels-normalized.tif");
                } else {
                    int indexCh = ArrayUtils.indexOf(channelNames, channels[0]);
                    imgVessels = BF.openImagePlus(options)[indexCh];
                }
                imgVessels.setLut(lut);
                IJ.run(imgVessels, "16-bit", "");
                
                ImagePlus imgMicro = null;
                if (!channels[1].equals("None")) {
                    int indexCh = ArrayUtils.indexOf(channelNames, channels[1]);
                    imgMicro = BF.openImagePlus(options)[indexCh];
                    imgMicro.setLut(lut);
                    IJ.run(imgMicro, "16-bit", "");
                }
                
                ImagePlus imgEndo = null;
                if (!channels[2].equals("None")) {                    
                    int indexCh = ArrayUtils.indexOf(channelNames, channels[2]);
                    imgEndo = BF.openImagePlus(options)[indexCh];
                    imgEndo.setLut(lut);
                    IJ.run(imgEndo, "16-bit", "");
                }
                
                // Load ROIs (if provided)
                tools.print("- Loading ROIs -");
                List<Roi> rois = tools.loadRois(imageDir + File.separator + rootName, imgVessels, rootName);
                
                // Segment vessels
                tools.print("- Segmenting vessels -");
                ImagePlus imgVesselMask = tools.vesselSegmentation(imgVessels, cal);
                
                // Compute vessels skeleton
                tools.print("- Skeletonizing vessels mask -");
                ImagePlus imgVesselSkel = tools.skeletonize3D(imgVesselMask, cal);
                // Prune vessels skeleton small branches
                imgVesselSkel = tools.pruneSkeleton(imgVesselSkel);
                
                // Compute vessels distance map
                tools.print("- Computing vessels distance map -");
                ImageFloat vesselDistMap = tools.distanceMap3D(imgVesselMask, false, cal);

                Objects3DIntPopulation microPop = new Objects3DIntPopulation();
                ImageFloat vesselDistMapInv = null;
                if (imgMicro != null) {
                    // Segment microglia
                    tools.print("- Segmenting microglia -");
                    microPop = tools.microSegmentation(imgMicro, cal);
                    // Compute vessels inverted distance map
                    tools.print("- Computing vessels inverted distance map -");
                    vesselDistMapInv = tools.distanceMap3D(imgVesselMask, true, cal);
                }
                
                Objects3DIntPopulation endoPop = new Objects3DIntPopulation();
                if (imgEndo != null) {
                    // Detect endothelial nuclei
                    tools.print("- Segmenting endothelial nuclei -");
                    endoPop = tools.endoSegmentation(imgEndo, cal);
                }
               
                // Save results
                tools.saveResults(rois, imgVesselMask, imgVesselSkel, vesselDistMap, vesselDistMapInv, microPop, endoPop, 
                                  imgVessels, imgMicro, imgEndo, cal, vesselResults, microResults, globalResults, rootName, outDir);

                tools.closeImage(imgVessels);
                tools.closeImage(imgVesselMask);
                tools.closeImage(imgVesselSkel);
                tools.closeImage(vesselDistMap.getImagePlus());
                if (imgMicro != null) tools.closeImage(imgMicro);
                if (imgMicro != null) tools.closeImage(vesselDistMapInv.getImagePlus());
                if (imgEndo != null) tools.closeImage(imgEndo);
            }
            globalResults.close();
            vesselResults.close();
            if (!channels[1].equals("None")) microResults.close();
        } catch (IOException ex) {
            Logger.getLogger(Vessels_Microglia_Endothelium.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(Vessels_Microglia_Endothelium.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("All done!");
    }
}
