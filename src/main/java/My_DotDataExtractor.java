/*
* Copyright (C) 2018 Federico Pevere and Benjamin Bruhn. All rights reserved.
*
* Permission to use, copy, modify, and distribute this software for any purpose without fee is hereby granted, provided that this entire notice is included in all copies of any software which is or includes a copy or modification of this software and in all copies of the supporting documentation for such software. Any for profit use of this software is expressly forbidden without first obtaining the explicit consent of the author.
*
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED WARRANTY. IN PARTICULAR, THE AUTHOR DOES NOT MAKE ANY REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
*
* Please put the following references in your publication if you use this plugin for your work:
*
* Text References:
*
* Pevere, F. , Bruhn, B. , Sangghaleh, F. , Hormozan, Y. , Sychugov, I. and Linnros, J. (2015), Effect of X‐ray irradiation on the blinking of single silicon nanocrystals. Phys. Status Solidi A, 212: 2692-2695.
*
* Bibtex:
*
* @article{blinkingSiQDs,
* author = {Pevere Federico and Bruhn Benjamin and Sangghaleh Fatemeh and Hormozan Yashar and Sychugov Ilya and Linnros Jan},
* title = {Effect of X‐ray irradiation on the blinking of single silicon nanocrystals},
* journal = {physica status solidi (a)},
* volume = {212},
* number = {12},
* pages = {2692-2695},
* doi = {10.1002/pssa.201532652}
* }
*
*/

//import packages (group of classes) that will be used in the plugin
//Java
import java.io.File;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Reader;
import java.io.PrintStream;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.*;
import java.awt.*;
import java.text.*;
//ImageJ
import ij.*;
import ij.gui.*;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.process.*;
import ij.plugin.frame.*;
import ij.plugin.filter.*;
import ij.plugin.frame.Editor;
import ij.io.*;
import ij.gui.Plot;
//statistics
//import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
//import org.apache.commons.math3.stat.StatUtils;

//plugin class

public class My_DotDataExtractor implements PlugInFilter {

    //BackgroundDot: private class which represents the background dot
    private class BackgroundDot {
	//class fields
        int bgX = 0; //x coordinate of the bg dot
        int bgY = 0; //y coordinate of the bg dot
        double bgDistance = 0.0; //bg dot distance
        double bgAngle = 0.0; //bg dot angle

        //class constructor
        BackgroundDot(int x, int y, MyQuantumDot referenceQD) {
            this.bgX = x; //assign x coord
            this.bgY = y; //assign y coord
            int dx = referenceQD.getX()-x; //calculate dx
            int dy = referenceQD.getY()-y; //calculate dy
            this.bgDistance = Math.sqrt(dx*dx+dy*dy); //assign distance from ref
            //angle calculation (in degrees)
            double angle = Math.abs(Math.atan(dy/dx) * 180 / Math.PI);
            // dx<0 && dy>0 top right from reference (NB x=y=0 corresponds to top left corner of the image)
            if(dx>0 && dy>0) { //top left from ref
                angle = 180 - angle;
            } else if (dx>0 && dy<0) { //bottom left from ref
                angle = 180 + angle;
            } else if (dx<0 && dy<0) { //bottom right from ref
                angle = 360 - angle;
            }
            this.bgAngle = angle; //assign angle calculated from ref
        }

        //class methods to get X and Y coordinates
        public int getX() {
            return this.bgX;
        }

        public int getY() {
            return this.bgY;
        }

        //class methods to get distance and angle calculated from referenceQD
        public double getDistance() {
            return this.bgDistance;
        }

        public double getAngle() {
            return this.bgAngle;
        }
    } //end of class BackgroundDot

    //BackgroundDotByDistance: private class which implements Comparator interface
    private class BackgroundDotsByDistance implements Comparator {
        //method  required to carry out contract with interface Comparator
        public int compare(Object obj1, Object obj2) {
            return (int) (((BackgroundDot)obj1).getDistance()-((BackgroundDot)obj2).getDistance());
        }
    } //end of class BackgroundDotsByDistance

    //BackgroundDotByDistanceAndAngle: private class (currently unused)
    private class BackgroundDotsByDistanceAndAngle implements Comparator {
        //method  required to carry out contract with interface Comparator
        @Override
        public int compare(Object obj1, Object obj2) {
            //how to take values from the calling function and previously calculated stuff into account here?!
            //only possible with global variables?
            //Modifications by Federico Pevere, 2015:
            //check the meaning of this function and its use in the code
            return 0;
        }
    } //end of class BackgroundDotByDistanceAndAngle

    ImagePlus  imp = null; //imp contains an ImageProcessor (2D image) or an ImageStack (3D)
    ImageStack stack = null; //stack represents an expandable array of images

    int nDots = 0; //initialize the number of the QDs
    int dotEdge = 5; //default edge length in px of the dot rectangle selection
    int bgDistance = 1; //bg px distance from edge of QD
    double frameTime = 1; //frame time in [s]
    boolean inputBgDotMap = false;
    int nBgDots = 0; //initialize the number of BG QDs

    int nHistogramBins = 50; //default number of bins of the integrated intensity histogram
    boolean plotHistogram = true; //default value

    int nBlinkingDots = 0; //initialize number of blinking dots

    Editor bugFinder = null; //create bugFinder editor

    //setup the plugin before use: an opened image stack is required otherwise
    //the plugin won't run (will display a msg in that case)

    /**
     *
     * @param arg
     * @param imp
     * @return
     */
    @Override
    public int setup(String arg, ImagePlus imp) {
        //IJ.register(My_DotDataExtractor.class);//probably not needed with
                                               //Java 1.2 or later
        this.imp = imp; //keep a reference of the image imp
        return DOES_ALL|STACK_REQUIRED; //vectors of binary flags that describes the plugin's properties
					//DOES_ALL: the filter handles all types of images
					//STACK_REQUIRED: the filter requires a stack
    }

    //run the plugin
    @Override
    public void run(ImageProcessor ip) {

        //information about plugin/debug
        bugFinder = new Editor(); //editor for code debugging
        bugFinder.display("BugFinder.txt","My_DotDataExtractor: an ImageJ plugin\n");
        bugFinder.display("BugFinder.txt","Authors: Federico Pevere (2015-2016, current developer) and Benjamin Bruhn (2012)\n");
        bugFinder.display("BugFinder.txt","Organization: KTH Royal Institute of Technology, Kista (Sweden)\n\n");
        bugFinder.display("BugFinder.txt","Contact (email): pevere@kth.se\n");
        bugFinder.display("BugFinder.txt","Acknowledgements: Swedish Research Council (VR), ADOPT\n\n");
        bugFinder.display("BugFinder.txt","Entering run method...\n\n");

        //image stack acquisition
        stack = imp.getStack(); //get the image stack (array of images)
        int stackWidth = ip.getWidth(); //get the width of the image stack
        int stackHeight = ip.getHeight(); //get the height of the image stack
        int stackSize = stack.getSize(); //get the size of the image stack
        int stackBitDepth =ip.getBitDepth(); //get the bit depth of the image stack

        //image stack information/debug
        bugFinder.append("Image stack info:\n");
        bugFinder.append("Size of the stack: "+stackSize+" images\n");
        bugFinder.append("Width of the image: "+stackWidth+" px\n");
        bugFinder.append("Height of the image: "+stackHeight+" px\n\n");
        bugFinder.append("Bitdepth of the image: "+stackBitDepth+" bit\n\n");

       //Future development:
       //do you want to 1. find the dots or 2. you have a previously saved dotmap?
       //1. Find the dots: calculate combined image (sum), run FindMaxima, set
       //appropriate threshold by asking approximate number of dots and then
       //save dotmap. 2. As usual.

        //read dotMap
        bugFinder.append("Loading dot map...\n");
        int[][] dotMap = readDotMap(); //load the dotMap
        bugFinder.append("...done!\n\n");
        bugFinder.append("Number of dots: "+nDots+"\n\n");

        //read bgDotMap
        int[][] bgDotMap = null;
        Frame parent = new Frame();
        YesNoCancelDialog dialogBgDotMap = new YesNoCancelDialog(parent,"Background dot map","Do you want to input a background dot map also?");
        if (!dialogBgDotMap.cancelPressed() && dialogBgDotMap.yesPressed()){
            //load the bgDotMap
            bugFinder.append("Loading bgDotMap[][]...\n");
            bgDotMap = readBgDotMap();
            bugFinder.append("...done!\n\n");
        }

        //input data extraction settings
        int gdDataExtrSettingsDigits=2;
        int gdDataExtrSettingsFieldWidth=5;
        bugFinder.append("Data extraction settings:\n\n");
        GenericDialog gdDataExtrSettings = new GenericDialog("Data extraction settings");
        gdDataExtrSettings.addNumericField("Dot selection edge",dotEdge,gdDataExtrSettingsDigits,gdDataExtrSettingsFieldWidth,"pixels");
        gdDataExtrSettings.addNumericField("Background distance from edge",bgDistance,gdDataExtrSettingsDigits,gdDataExtrSettingsFieldWidth,"pixels");
        gdDataExtrSettings.showDialog();
        //check validity of input
        if (gdDataExtrSettings.wasCanceled()){
            IJ.showMessage("Default values will be used:\n "
                          + "Dot selection edge: "+dotEdge+" px\n"
                          + "Background distance from edge: "+bgDistance+" px");
        }else{
            //selection edge
            int tempDotEdge = (int) gdDataExtrSettings.getNextNumber();
            if (!gdDataExtrSettings.invalidNumber() && tempDotEdge>0){
                dotEdge = tempDotEdge;
            }else{
                IJ.showMessage("Invalid number inserted! The selection area of the QD will be of "+dotEdge+" px edge. \n\n");
            }
            bugFinder.append("Edge: "+dotEdge+" px\n");
            //bg distance
            int tempBgDistance = (int) gdDataExtrSettings.getNextNumber();
            if (!gdDataExtrSettings.invalidNumber() && tempBgDistance>0){
                bgDistance = tempBgDistance;
            }else{
                IJ.showMessage("Invalid numbers inserted! The background will be evaluated "+bgDistance+" px far from the edge of the QD selection area. \n\n");
            }
            bugFinder.append("BG distance: "+bgDistance+" px\n\n");
        }

        //input frame time
        int gdFrameTimeDigits=2;
        int gdFrameTimeWidth=5;
        bugFinder.append("Experimental acquisition settings:\n\n");
        GenericDialog gdFrameTime = new GenericDialog("Experimental acquisition settings");
        gdFrameTime.addNumericField("Frame time",frameTime,gdFrameTimeDigits,gdFrameTimeWidth,"s");
        gdFrameTime.showDialog();
        //check validity of input
        if (gdFrameTime.wasCanceled()){
            IJ.showMessage("Default value will be used:\n "
                          + "Frame time: "+frameTime+" s");
        }else{
            double tempFrameTime = gdFrameTime.getNextNumber();
            if (!gdFrameTime.invalidNumber() && tempFrameTime>0){
                frameTime = tempFrameTime;
            }else{
                IJ.showMessage("Invalid number inserted! The frame time will be "+frameTime+" s. \n\n");
            }
        }
        bugFinder.append("Frame time: "+frameTime+" s\n");

        //create dot list (collection class containing the dots) with an initial
        //capacity of 10
        ArrayList<MyQuantumDot> dotList = new ArrayList<MyQuantumDot>();
        //extract the integrated intensities
        bugFinder.append("Extracting QDs intensity traces:\n\n");
        IJ.showStatus("Extracting intensity traces..");
        IJ.showProgress(0.0); //updates the progress bar, where 0<=progress<=1.0
        //fill in the dot list
        for(int i=0; i<nDots; i++) {
            //add new QD
            //bugFinder.append("Dot "+(i+1)+" @ ("+dotMap[i][0]+","+dotMap[i][1]+") px..");
            dotList.add(new MyQuantumDot(dotMap[i][0], dotMap[i][1], dotEdge, bgDistance, imp));
            IJ.showProgress((double) (i+1)/nDots); //update status bar
            //bugFinder.append("..intensity trace extracted!\n");
            //bugFinder.append("Integrated intensity in 1st frame: "+dotList.get(i).getIntegratedIntensity(0)+" arb. u.\n\n");
            //calculate the background intensity when the bgDotMap is not empty
            if(bgDotMap!=null) {
                setBgIntegrIntens(dotList.get(i), bgDotMap);
            }
        }

        bugFinder.append("..done!\n");

        //mandatory output:
        //display dotMap
        bugFinder.append("Displaying the dot map...\n");
        displayOverlay(dotList);
        bugFinder.append("...done!\n\n");
        //save all the integrated intensities
        bugFinder.append("Saving all integrated intensities to file...\n");
        saveAllIntegrIntensToFile(dotList,stackSize);
        bugFinder.append("...done!\n\n");

        //ASK WHAT TO PLOT

        //create time array
        double[][] time = new double[nDots][stackSize];
        for(int dotNo=0; dotNo<nDots; dotNo++){
            for(int frameCounter=0; frameCounter<stackSize;frameCounter++){
                time[dotNo][frameCounter]=(frameCounter+1)*frameTime;
            }
        }

        //create integrated intensities arrays
        double[][] integrIntens = new double[nDots][stackSize];
        double[][] bgIntegrIntens = new double[nDots][stackSize];
        double[][] bgCorrIntegrIntens = new double[nDots][stackSize];

        for(int dotNo=0; dotNo < nDots; dotNo++){
            MyQuantumDot qd = dotList.get(dotNo);
            integrIntens[dotNo]= qd.getIntegrIntens();
            bgIntegrIntens[dotNo]= qd.getBgIntegrIntens();
            bgCorrIntegrIntens[dotNo]= qd.getBgCorrIntegrIntens();
        }

        //Ask what integr. intens. to plot
        String[] gdPlotChooseIntegrIntensLabels = {"Integr. Intens. Signal+BG","Integr. Intens. BG","Integr. Intens. Signal"};
        boolean[] gdPlotChooseIntegrIntensDefaultValues = {true, true, true};
        int gdPlotChooseIntegrIntensColumns=1;
        int gdPlotChooseIntegrIntensDigits=1;
        int gdPlotChooseIntegrIntensFieldWidth=5;
        int gdPlotChooseIntegrIntensFrameWidth=1000;
        int gdPlotChooseIntegrIntensFrameHeight=500;
        GenericDialog gdPlotChooseIntegrIntens = new GenericDialog("Plot details");
        gdPlotChooseIntegrIntens.addCheckboxGroup(gdPlotChooseIntegrIntensLabels.length, gdPlotChooseIntegrIntensColumns, gdPlotChooseIntegrIntensLabels, gdPlotChooseIntegrIntensDefaultValues);
        gdPlotChooseIntegrIntens.addNumericField("Plots width", gdPlotChooseIntegrIntensFrameWidth, gdPlotChooseIntegrIntensDigits, gdPlotChooseIntegrIntensFieldWidth, "pixels");
        gdPlotChooseIntegrIntens.addNumericField("Plots height", gdPlotChooseIntegrIntensFrameHeight, gdPlotChooseIntegrIntensDigits, gdPlotChooseIntegrIntensFieldWidth, "pixels");
        gdPlotChooseIntegrIntens.showDialog();

        if (!gdPlotChooseIntegrIntens.wasCanceled() && gdPlotChooseIntegrIntens.wasOKed()){

            int tempPlotWidth = (int) gdPlotChooseIntegrIntens.getNextNumber();
            int tempPlotHeight = (int) gdPlotChooseIntegrIntens.getNextNumber();

            if(!gdPlotChooseIntegrIntens.invalidNumber() &&  tempPlotWidth>0 && tempPlotHeight>0){
                gdPlotChooseIntegrIntensFrameWidth=tempPlotWidth;
                gdPlotChooseIntegrIntensFrameHeight=tempPlotHeight;
            }

            //bugFinder.append("Plotting all the integrated intensities...\n");
            if(gdPlotChooseIntegrIntens.getNextBoolean()){
            //plot integrated intensity traces
            bugFinder.append("Creating plot stack for integr. intens. traces...\n");
            ImagePlus plotStackIntegrIntensImp = getPlotStack(time,integrIntens,"Integrated intensities (signal+background)","QD","Time [s]","Integr. intens. [arb. u.]",gdPlotChooseIntegrIntensFrameWidth,gdPlotChooseIntegrIntensFrameHeight);
            plotStackIntegrIntensImp.show();
            bugFinder.append("...done!\n\n");
            }

            if(gdPlotChooseIntegrIntens.getNextBoolean()){
            //plot background integrated intensity traces
            bugFinder.append("Creating plot stack for bg integr. intens. traces...");
            ImagePlus plotStackBgIntegrIntensImp = getPlotStack(time,bgIntegrIntens,"Background "+"(background at "+bgDistance+" px far from dot edge)","QD","Time [s]","Integr. intens. [arb. u.]",gdPlotChooseIntegrIntensFrameWidth,gdPlotChooseIntegrIntensFrameHeight);
            plotStackBgIntegrIntensImp.show();
            bugFinder.append("...done!\n\n");
            }

            if(gdPlotChooseIntegrIntens.getNextBoolean()){
            //plot integrated intensity traces
            bugFinder.append("Creating plot stack for bg-corr. integr. intens. traces...");
            ImagePlus plotStackBgCorrIntegrIntensImp = getPlotStack(time,bgCorrIntegrIntens,"Integrated intensities (signal)","QD","Time [s]","Integr. intens. [arb. u.]",gdPlotChooseIntegrIntensFrameWidth,gdPlotChooseIntegrIntensFrameHeight);
            plotStackBgCorrIntegrIntensImp.show();
            bugFinder.append("...done!\n\n");
            }

//            // uncomplete code
//            bugFinder.append("Entering test plot...\n");
//            // Set up stack
//            int dotNo=0;
//            // Plot first plot just to get dimensions
//            double stddevY=dotList.get(dotNo).getStddevBgCorrIntegrIntens();
//            double meanY=dotList.get(dotNo).getMeanBgCorrIntegrIntens();
//
//            double minY=(meanY-4*stddevY) > dotList.get(dotNo).getMinBgCorrIntegrIntens()?(meanY-4*stddevY):dotList.get(dotNo).getMinBgCorrIntegrIntens();
//            double maxY=(meanY+4*stddevY) < dotList.get(dotNo).getMaxBgCorrIntegrIntens()?(meanY+4*stddevY):dotList.get(dotNo).getMaxBgCorrIntegrIntens();
//
//            Plot tempPlot = new Plot ("Test", "Time [s]", "Integr. Int.", time[dotNo], bgCorrIntegrIntens[dotNo]);
//            tempPlot.setLimits(0, StatUtils.max(time[dotNo]), minY, maxY);
//            ImageProcessor tempIp = tempPlot.getProcessor();
//            int plotStackWidth = tempIp.getWidth();
//            bugFinder.append("stackWidth: "+stackWidth+"\n");
//            int plotStackHeight = tempIp.getHeight();
//            bugFinder.append("stackHeigth: "+stackHeight+"\n");
//
//            // Set up stack
//            ImageStack plotStack = new ImageStack(plotStackWidth, plotStackHeight);
//
//            bugFinder.append("ImageStack created...\n");
//
//            // Add each plot
//            bugFinder.append("Entering for loop...\n");
//            IJ.showProgress(0.0);
//            for (dotNo=0; dotNo<nDots; dotNo++) {
//                IJ.showStatus("Creating Test plot stack...");
//                // Plot intensity graph
//                stddevY=Math.sqrt(StatUtils.variance(bgCorrIntegrIntens[dotNo]));
//                meanY=StatUtils.mean(bgCorrIntegrIntens[dotNo]);
//                minY=(meanY-4*stddevY) > StatUtils.min(bgCorrIntegrIntens[dotNo])?(meanY-4*stddevY):StatUtils.min(bgCorrIntegrIntens[dotNo]);
//                maxY=(meanY+4*stddevY) < StatUtils.max(bgCorrIntegrIntens[dotNo])?(meanY+4*stddevY):StatUtils.max(bgCorrIntegrIntens[dotNo]);
//                Plot plot = new Plot ("Test", "Time [s]", "Integr. Intens. [arb. u.]", time[dotNo], bgCorrIntegrIntens[dotNo]);
//                plot.setLimits(0, StatUtils.max(time[dotNo]), minY, maxY);
//                //plot.setLineWidth(1);
//                //plot.setColor(Color.red);
//                //plot.changeFont(new Font("Helvetica", Font.PLAIN, 16));
//                //plot.setFrameSize(plotStackWidth, plotStackHeight);
//                //plot.show();
//                bugFinder.append("Plot created...\n");
//                ImageProcessor ipPlot = plot.getProcessor();
//                bugFinder.append("ipPlot created...\n");
//                plotStack.addSlice("Test plot of dot "+(dotNo+1), ipPlot);
//                bugFinder.append("Test plot of Dot"+(dotNo+1)+"plotted...\n");
//                IJ.showProgress((double) (dotNo+1)/bgCorrIntegrIntens.length);
//            }
//
//        ImagePlus testPlot = new ImagePlus ("test plot stack", plotStack);
//        testPlot.show();
//        bugFinder.append("End of test plot...\n");
        }

        //do you want to save  integr. intens.?
        //YesNoCancelDialog dialogStatAllIntegrIntens = new YesNoCancelDialog(parent,"Statistics of all integr. intens.","Save the statistics of all (signal, bg,...) the integr. intens.?");
        //if (!dialogStatAllIntegrIntens.cancelPressed() && dialogStatAllIntegrIntens.yesPressed()){

        //save additional info
        //bugFinder.append("Saving statistics of all integr. intens. to file...\n");
        //saveStatAllIntegrIntensToFile(dotList);
        //bugFinder.append("...done!\n\n");
        //save histograms
//        bugFinder.append("Saving integrated intensity histograms to file...\n");
//        saveIntegratedIntensityHistogramsToFile(dotList,nHistogramBins);
//        bugFinder.append("...done!\n\n");
        //}

//
//show the results of the fitting and of the blinking selection algorithm
//
//        bugFinder.append("Results of the blinking detection algorithm:\n\n");
//        DecimalFormat df = new DecimalFormat("####0.00");
//
//        for(int i=0; i< nDots; i++){
//            //debug
//            bugFinder.append("Dot"+dotList.get(i).ID+":\n");
//            bugFinder.append("Initial fitting parameters:\n");
//            bugFinder.append("a (areaON*sqrt(2*pi)*sdON): "+dotList.get(i).getInitialFittingParam()[0]+"\n");
//            bugFinder.append("b (meanONvalue): "+dotList.get(i).getInitialFittingParam()[1]+"\n");
//            bugFinder.append("c (sdON): "+dotList.get(i).getInitialFittingParam()[2]+"\n");
//            bugFinder.append("d (areaOFF*sqrt(2*pi)*sdOFF): "+dotList.get(i).getInitialFittingParam()[3]+"\n");
//            bugFinder.append("e (meanOFFvalue): "+dotList.get(i).getInitialFittingParam()[4]+"\n");
//            bugFinder.append("f (sdOFF): "+dotList.get(i).getInitialFittingParam()[5]+"\n");
//            bugFinder.append("Fit goodness: "+dotList.get(i).cf.getFitGoodness()+"\n");
//            bugFinder.append(dotList.get(i).cf.getResultString()+"\n\n");
//            if(dotList.get(i).getTwoLevelBlinking()){
//                nBlinkingDots++;
//                bugFinder.append("Dot"+(i+1)+" is a BLINKING DOT !!!\n");
//            } else
//                bugFinder.append("Dot"+(i+1)+" is NOT a BLINKING DOT !!!\n");
//
//                bugFinder.append("Integrated intensity mean ON level: "+df.format(dotList.get(i).getIntIntnON())+" arb. u.\n");
//                bugFinder.append("Integrated intensity sd ON level: "+df.format(dotList.get(i).getSdIntIntnON())+" arb. u.\n");
//                bugFinder.append("Integrated intensity mean OFF level: "+df.format(dotList.get(i).getIntIntnOFF())+" arb. u.\n");
//                bugFinder.append("Integrated intensity sd OFF level: "+df.format(dotList.get(i).getSdIntIntnOFF())+" arb. u.\n");
//                bugFinder.append("Integrated intensity mean amplitude (ONtoOFF): "+df.format(dotList.get(i).getIntIntnONtoOFF())+" arb. u.\n");
//                bugFinder.append("Integrated intensity sd amplitude (ONtoOFF) : "+df.format(dotList.get(i).getSdIntIntnONtoOFF())+" arb. u.\n");
//                bugFinder.append("Duty cycle (ON): "+df.format(dotList.get(i).getDutyCyleFromHistogram()*100)+"% \n");
//                bugFinder.append("Automatic ON-OFF threshold: "+df.format(dotList.get(i).getOnOffThreshold())+" arb. u.\n");
//                bugFinder.append("\n");
//
//        }
//
//        //analysis of blinking QDs
//        IJ.showMessage("Blinking search result:\n","Number of BLINKING dots: "+nBlinkingDots+"\n\n");
//        bugFinder.append("In total, "+nBlinkingDots+" BLINKING dots were found. Namely:\n");
//        for(int i=0; i<nDots; i++){
//            if(dotList.get(i).getTwoLevelBlinking()){
//                bugFinder.append("Dot "+dotList.get(i).getID()+"\n");
//            }
//        }
//        bugFinder.append("\n");
//
//        //Analysis of BLINKING dots
//        bugFinder.append("Analysis of BLINKING dots...(work in progress)\n\n");
//        Frame parent1 = new Frame();
//        YesNoCancelDialog dialogBlinkingAnalysis = new YesNoCancelDialog(parent1,"Blinking analysis","Do you want to analyze BLINKING dots?");
//
//        if (!dialogBlinkingAnalysis.cancelPressed() && dialogBlinkingAnalysis.yesPressed()){
//            //get ON and OFF times, blinking frequency, ON- and OFF- time constants
//            //save results
//
//            //saveOnOffTimes();
//            //saveBlinkingParameters();
//
//        }
//
//
        bugFinder.append("Leaving run.\n\n");
    }

    public void setBgIntegrIntens(MyQuantumDot qd, int[][] bgDotMap) {

        //determine suitable background sources
        ArrayList<BackgroundDot> bgDotList = new ArrayList<BackgroundDot>();
        for(int i=0; i<nBgDots; i++) {
            bgDotList.add(new BackgroundDot(bgDotMap[i][0],bgDotMap[i][1],qd));
        }
        Collections.sort(bgDotList,new BackgroundDotsByDistance());
        BackgroundDot backgroundDot1 = bgDotList.get(0);

        //for each frame calculate background and set in QD
        int sum = 0;
        int[] pixels = new int[dotEdge*dotEdge];
        ImageProcessor ip = null;
        for(int i=0; i<stack.getSize(); i++) {
            ip = stack.getProcessor(i+1);
            for(int j=0; j<dotEdge; j++){
                for(int k=0; k<dotEdge; k++){
                    pixels[k+dotEdge*j] = ip.getPixel(backgroundDot1.getX()-(int)((dotEdge-1)/2)+k,backgroundDot1.getY()-(int)((dotEdge-1)/2)+j);
                }
            }
            for(int h=0; h<dotEdge*dotEdge; h++){
                sum += pixels[h];
            }
            qd.setBgIntegrIntens(i,sum);
            sum=0;
        }
    }

    ImagePlus getPlotStack (double[][] xData, double[][] yData, String stackLabel, String plotTitle, String xLabel, String yLabel, int width, int height){
        //bugFinder.append("Entering getPlotStack...\n");
        int plotNo=0;
        // Plot first plot just to get dimensions

        Plot tempPlot = new Plot (plotTitle, xLabel, yLabel, xData[plotNo], yData[plotNo]);
        ImageProcessor tempIp = tempPlot.getProcessor();
        int stackWidth = tempIp.getWidth();
        //bugFinder.append("stackWidth: "+stackWidth+"\n");
        int stackHeight = tempIp.getHeight();
        //bugFinder.append("stackHeigth: "+stackHeight+"\n");

        // Set up stack
        ImageStack plotStack = new ImageStack(stackWidth, stackHeight);

        // Add each plot
        IJ.showProgress(0.0);
	for (plotNo=0; plotNo<yData.length; plotNo++) {
            IJ.showStatus("Creating "+stackLabel+" plot stack...");
            // Plot intensity graph
            Plot plot = new Plot (plotTitle, xLabel, yLabel, xData[plotNo], yData[plotNo]);
            //plot.setSize(width, height); //set canvas size
            ImageProcessor ip = plot.getProcessor();
            plotStack.addSlice(plotTitle+" "+(plotNo+1), ip);
            IJ.showProgress((double) (plotNo+1)/nDots);
	}

        return new ImagePlus (stackLabel, plotStack);
    }

//    ImagePlus getPlotIntIntensStack (int stackSize, double frameTime, ArrayList<MyQuantumDot> dotList) {
//
//        double[] xAxis = new double[stackSize];
//        for(int i=0; i<stackSize; i++) {
//            xAxis[i] = (double) (i+1)*frameTime;
//        }

//        double[] yAxis = Arrays.stream(dotList.get(0).getIntegratedIntensity()).asDoubleStream().toArray();
//
//        Plot plot = new Plot ("Trace for dot " + dotID, "frame no.", "Int. intensity (a.u.)",xAxis,yAxisD);
//
//        // Plot first quantum dot to get dimensions
//        ImageProcessor tempIp = getPlotIntIntensObject(xAxis,dotList.get(0).getIntegratedIntensity(),0,stackSize).getProcessor();
//        int stackWidth = tempIp.getWidth();
//        int stackHeight = tempIp.getHeight();
//
//        // Set up stack
//        ImageStack plotStack = new ImageStack(stackWidth, stackHeight);
//
//        // Add each plot
//        IJ.showProgress(0.0);
//	for (int dotNo=0; dotNo<nDots; dotNo++) {
//            IJ.showStatus("Creating int. intens. plot stack...");
//            ImageProcessor ip = getPlotIntIntensObject(xAxis,dotList.get(dotNo).getIntegrIntens(),(dotNo+1)).getProcessor();
//            plotStack.addSlice("dot " + (dotNo+1), ip);
//            IJ.showProgress((double) (dotNo+1)/nDots);
//	}
//
//	return new ImagePlus ("QD Integrated intensities", plotStack);
//    }
//
//    Plot getPlotIntIntensObject (double[] xAxis, int[] yAxis, int dotID) {
//
//        double[] yAxisD = Arrays.stream(yAxis).asDoubleStream().toArray();
//        // Plot intensity graph
//        Plot plot = new Plot ("Trace for dot " + dotID, "frame no.", "Int. intensity (a.u.)",xAxis,yAxisD);
//        return plot;
//    }

//    ImagePlus getPlotBgIntIntensStack (int stackSize, ArrayList<MyQuantumDot> dotList) {
//
//        double[] xAxis = new double[stackSize];
//        for(int i=0; i<stackSize; i++) {
//            xAxis[i] = (double) (i+1);
//        }
//
//        // Plot first quantum dot to get dimensions
//        ImageProcessor tempIp = getPlotBgIntIntensObject(xAxis,dotList.get(0).getBackgroundIntensity(),0,stackSize).getProcessor();
//        int stackWidth = tempIp.getWidth();
//        int stackHeight = tempIp.getHeight();
//
//        // Set up stack
//        ImageStack plotStack = new ImageStack(stackWidth, stackHeight);
//
//        // Add each plot
//        IJ.showProgress(0.0);
//	for (int dotNo=0; dotNo<nDots; dotNo++) {
//            IJ.showStatus("Creating background int. intens. plot stack...");
//            ImageProcessor ip = getPlotBgIntIntensObject(xAxis,dotList.get(dotNo).getBackgroundIntensity(),dotList.get(dotNo).getID(),stackSize).getProcessor();
//            plotStack.addSlice("dot " + dotList.get(dotNo).getID(), ip);
//            IJ.showProgress((double) (dotNo+1)/nDots);
//	}
//
//	return new ImagePlus ("QD Bacground intensities", plotStack);
//    }
//
//    Plot getPlotBgIntIntensObject (double[] xAxis, int[] yAxis, int dotID, int stackSize) {
//
//    double[] yAxisD = new double[stackSize];
//    for(int i=0; i<stackSize; i++) {
//        yAxisD[i] = (double) yAxis[i];
//    }
//
//    // Plot intensity graph
//    Plot plot = new Plot ("Trace for dot " + dotID, "frame no.", "Bg int. intensity (a.u.)",xAxis,yAxisD);
//    return plot;
//    }
//
//    ImagePlus getPlotBgCorrIntIntensStack (int stackSize, ArrayList<MyQuantumDot> dotList) {
//
//        double[] xAxis = new double[stackSize];
//        for(int i=0; i<stackSize; i++) {
//            xAxis[i] = (double) (i+1);
//        }
//
//        // Plot first quantum dot to get dimensions
//        ImageProcessor tempIp = getPlotBgCorrIntIntensObject(xAxis,dotList.get(0).getBgCorrectedIntIntensity(),0,stackSize).getProcessor();
//        int stackWidth = tempIp.getWidth();
//        int stackHeight = tempIp.getHeight();
//
//        // Set up stack
//        ImageStack plotStack = new ImageStack(stackWidth, stackHeight);
//
//        // Add each plot
//        IJ.showProgress(0.0);
//	for (int dotNo=0; dotNo<nDots; dotNo++) {
//            IJ.showStatus("Creating bg-corrected int. intens. plot stack...");
//            ImageProcessor ip = getPlotBgCorrIntIntensObject(xAxis,dotList.get(dotNo).getBgCorrectedIntIntensity(),dotList.get(dotNo).getID(),stackSize).getProcessor();
//            plotStack.addSlice("dot " + dotList.get(dotNo).getID(), ip);
//            IJ.showProgress((double) (dotNo+1)/nDots);
//	}
//
//	return new ImagePlus ("QD Integrated intensities (bg-corrected)", plotStack);
//    }
//
//    Plot getPlotBgCorrIntIntensObject (double[] xAxis, int[] yAxis, int dotID, int stackSize) {
//
//    double[] yAxisD = new double[stackSize];
//    for(int i=0; i<stackSize; i++) {
//        yAxisD[i] = (double) yAxis[i];
//    }
//
//    // Plot intensity graph
//    Plot plot = new Plot ("Trace for dot " + dotID, "frame no.", "Int. intensity (bg subtracted) (a.u.)",xAxis,yAxisD);
//    return plot;
//    }
//
//
//    ImagePlus getHistogramStack (ArrayList<MyQuantumDot> dotList) {
//
//        // Plot first quantum dot to get dimensions
//        ImageProcessor tempIp = getHistogramObject(dotList.get(0),0).getProcessor();
//        int stackWidth = tempIp.getWidth();
//        int stackHeight = tempIp.getHeight();
//
//        // Set up stack
//        ImageStack plotStack = new ImageStack(stackWidth, stackHeight);
//
//        // Add each plot
//        IJ.showProgress(0.0);
//                for (int dotNo=0; dotNo<nDots; dotNo++) {
//                    IJ.showStatus("Creating histogram stack...");
//                    ImageProcessor ip = getHistogramObject(dotList.get(dotNo),dotList.get(dotNo).getID()).getProcessor();
//                    plotStack.addSlice("dot " + dotList.get(dotNo).getID(), ip);
//                    IJ.showProgress((double) (dotNo+1)/nDots);
//                }
//
//            return new ImagePlus ("Quantum dot histograms", plotStack);
//    }
//
//
//    Plot getHistogramObject (MyQuantumDot dot, int dotID) {
//
//        double histogramRange = dot.getIntegratedIntensityHistogram().getMaxBinEnd()-dot.getIntegratedIntensityHistogram().getMinBinStart();
//        int Nbins = dot.getIntegratedIntensityHistogram().getFreq().length;
//        double histogramBinSize = histogramRange/Nbins;
//
//        double[] xAxis = new double[Nbins];
//        double[] yAxis = new double[Nbins];
////        double[] xAxisFit = dot.cf.getXPoints();
////        double[] yAxisFit = new double[Nbins];
//
//        for(int i = 0; i < Nbins; i++) {
//            xAxis[i] = (double) dot.getIntegratedIntensityHistogram().getFreq()[i];
//            yAxis[i] =  dot.getIntegratedIntensityHistogram().getMinBinStart()+histogramBinSize/2+i*histogramBinSize;
//        }
//
//        // Plot intensity graph
//        //Plot plot = new Plot ("Histogram of dot " + dotID, "Occurence", "Intensity (arb. u.)", xAxis, yAxis);
//        //plot.setLineWidth (1);
//
//        Plot plot = new Plot ("Histogram of dot " + dotID, "Occurence", "Intensity (arb. u.)",xAxis,yAxis);
//        //int symbol = Plot.CIRCLE;
//	//plot.addPoints (xAxis, yAxisD, symbol);
//        //Color red = new Color(1,0,0);
//        //plot.setColor(red);
//        //plot.setLineWidth (1);
//        //int symbol = Plot.CIRCLE;
//        //plot.addPoints (xAxis, yAxis, symbol);
//
//        return plot;
//    }

    void saveAllIntegrIntensToFile(ArrayList<MyQuantumDot> dotList, int stackSize) {
        //save integr. intens. of the dots
        SaveDialog od = new SaveDialog("Select the file to save the integr. intens. data in", "qdIntegrIntens",".dat");
        if(od.getFileName()==null) {
            return;
        }
        String inputPath = od.getDirectory();
        String chosenFile = od.getFileName();
        File outFile = new File(inputPath + chosenFile);
        //save integr. intens. of the bg dots
        SaveDialog odBg = new SaveDialog("Select the file to save the bg integr. intens. data in", "bgIntegrIntens",".dat");
        if(odBg.getFileName()==null) {
            return;
        }
        String inputPathBG = odBg.getDirectory();
        String chosenFileBG = odBg.getFileName();
        File outFileBG = new File(inputPathBG + chosenFileBG);
        //save bg-corrected integr. intens. of the dots
        SaveDialog odBgCorrIntegrIntens = new SaveDialog("Select the file to save the bg-corrected integr. intens. data in", "qdBgCorrIntegrIntens",".dat");
        if(odBgCorrIntegrIntens.getFileName()==null) {
            return;
        }
        String inputPathBgCorrIntegrIntens = odBgCorrIntegrIntens.getDirectory();
        String chosenFileBgCorrIntegrIntens = odBgCorrIntegrIntens.getFileName();
        File outFileBgCorrIntegrIntens = new File(inputPathBgCorrIntegrIntens + chosenFileBgCorrIntegrIntens);
        //bugFinder.append("opened file\n");
        try {
            PrintStream psIntegrIntens = new PrintStream(outFile);
            PrintStream psBgIntegrIntens = new PrintStream(outFileBG);
            PrintStream psBgCorrIntegrIntens = new PrintStream(outFileBgCorrIntegrIntens);
            //bugFinder.append("opened printstream\n");
            IJ.showProgress(0.0);
            IJ.showStatus("Saving extracted data to file...");
            //insert "Frame" column label
            psIntegrIntens.print("Frame \t ");
            psBgIntegrIntens.print("Frame \t ");
            psBgCorrIntegrIntens.print("Frame \t ");
            //insert "Time" column label
            psIntegrIntens.print("Time \t ");
            psBgIntegrIntens.print("Time \t ");
            psBgCorrIntegrIntens.print("Time \t ");
            //insert column labels: QD 1 QD 2 ... (dot number)
            for(int dotNo=0; dotNo<nDots; dotNo++){
                psIntegrIntens.print("QD "+(dotNo+1)+"\t");
                psBgIntegrIntens.print("QD "+(dotNo+1)+"\t");
                psBgCorrIntegrIntens.print("QD "+(dotNo+1)+"\t");
            }
            //new line and skip first two columns
            psIntegrIntens.print("\n\t\t");
            psBgIntegrIntens.print("\n\t\t");
            psBgCorrIntegrIntens.print("\n\t\t");
            //insert column labels: (100,360) px ... (dot coordinates)
            for(int dotNo=0; dotNo<nDots; dotNo++){
                psIntegrIntens.print("("+dotList.get(dotNo).centerX+","+dotList.get(dotNo).centerY+") px\t"); //dot coordinates
                psBgIntegrIntens.print("("+dotList.get(dotNo).centerX+","+dotList.get(dotNo).centerY+") px\t");
                psBgCorrIntegrIntens.print("("+dotList.get(dotNo).centerX+","+dotList.get(dotNo).centerY+") px\t");
            }
            //new line
            psIntegrIntens.print("\n");
            psBgIntegrIntens.print("\n");
            psBgCorrIntegrIntens.print("\n");
            //add units
            psIntegrIntens.print("frames\t");
            psBgIntegrIntens.print("frames\t");
            psBgCorrIntegrIntens.print("frames\t");
            psIntegrIntens.print("s\t");
            psBgIntegrIntens.print("s\t");
            psBgCorrIntegrIntens.print("s\t");
            //insert column labels: (100,360) px ... (dot coordinates)
            for(int dotNo=0; dotNo<nDots; dotNo++){
                psIntegrIntens.print("arb. u.\t"); //dot coordinates
                psBgIntegrIntens.print("arb. u.\t");
                psBgCorrIntegrIntens.print("arb. u.\t");
            }
            //new line and skip first two columns
            psIntegrIntens.print("\n");
            psBgIntegrIntens.print("\n");
            psBgCorrIntegrIntens.print("\n");

            double time =0;

            for(int frame=0; frame<stackSize; frame++) {
                if(IJ.escapePressed()) {
                    psIntegrIntens.close();
                    psBgIntegrIntens.close();
                    psBgCorrIntegrIntens.close();
                    return;
                }
                //print number of frame
                psIntegrIntens.print(frame+1+"\t");
                psBgIntegrIntens.print(frame+1+"\t");
                psBgCorrIntegrIntens.print(frame+1+"\t");
                //print time in [s]
                time = (frame+1)*frameTime;
                psIntegrIntens.print(time+"\t");
                psBgIntegrIntens.print(time+"\t");
                psBgCorrIntegrIntens.print(time+"\t");

                //bugFinder.append("Saving info from frame "+frame+"\n");
                for(int dotNo=0; dotNo<nDots; dotNo++){
                    //IJ.showStatus("Saving extracted data to file...");
                    psIntegrIntens.print(dotList.get(dotNo).getIntegrIntens(frame));
                    //bugFinder.append("wrote intensity of dot "+dotNo+"\n");
                    psIntegrIntens.print("\t");
                    psBgIntegrIntens.print(dotList.get(dotNo).getBgIntegrIntens(frame));
                    //bugFinder.append("wrote intensity of dot "+dotNo+"\n");
                    psBgIntegrIntens.print("\t");
                    //bugFinder.append("wrote tab\n");
                    //IJ.showProgress((double) (frame+1)/stackSize);
                    psBgCorrIntegrIntens.print(dotList.get(dotNo).getIntegrIntens(frame)-dotList.get(dotNo).getBgIntegrIntens(frame));
                    psBgCorrIntegrIntens.print("\t");
                }
                IJ.showProgress((double) (frame+1)/stackSize);
                //new line
                psIntegrIntens.print("\n");
                psBgIntegrIntens.print("\n");
                psBgCorrIntegrIntens.print("\n");
            }
            psIntegrIntens.close();
            psBgIntegrIntens.close();
            psBgCorrIntegrIntens.close();
        } catch(Exception e) {
            bugFinder.append("Exception in opening and saving files\n");
        }
    }

//    // uncomplete code
//    void saveStatAllIntegrIntensToFile(ArrayList<MyQuantumDot> dotList) {
//        SaveDialog od = new SaveDialog("Select the file to save the statistic of all integr. intens. in", "statAllIntegrIntens",".dat");
//        if(od.getFileName()==null) {
//            return;
//        }
//        String inputPath = od.getDirectory();
//        String chosenFile = od.getFileName();
//        File outFile = new File(inputPath + chosenFile); //new File object call outFile
//        bugFinder.append("opened file\n");
//        try {
//            IJ.showProgress(0.0);
//            PrintStream psStatAllIntegrIntens = new PrintStream(outFile);
//            bugFinder.append("opened printstream\n");
//            IJ.showProgress(0.0);
//
//            //header
//            psStatAllIntegrIntens.println("Dot \t Mean Integr. Intens. \t Stddev Integr. Intens. \t"
//                    + "Max Integr. Intens. \t Min Integr. Intens. \t Mean Bg Integr. Intens. \t"
//                    + "Stddev Bg Integr. Intens. \t Max Bg Integr. Intens. \t Min Bg Integr. Intens \t"
//                    + "Mean Bg-corr Integr. Intens. \t Stddev Bg-corr Integr. Intens. \t"
//                    + "Max Bg-corr Integr. Intens. \t Min Bg-corr Integr. Intens.\n");
//
//            //data
//            for(int dotNo=0; dotNo<nDots; dotNo++) {
//                IJ.showProgress((double) (dotNo+1)/nDots);
//                if(IJ.escapePressed()) {
//                    psStatAllIntegrIntens.close();
//                    return;
//                }
//                bugFinder.append("Saving info from dot "+dotNo+"\n");
//                MyQuantumDot qd = dotList.get(dotNo);
//                psIntegrIntens.print(qd.getID());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getMaxIntegratedIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getMinIntegratedIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getMaxMinDifferenceIntegratedIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getMeanIntegratedIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getStandardDeviationIntegratedIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getMeanBackgroundIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getSdBackgroundIntensity());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getCurveFitter().getFormula());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getCurveFitter().getFitGoodness());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getIntIntnON());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getSdIntIntnON());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getIntIntnOFF());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getSdIntIntnOFF());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getTwoLevelBlinking()?"Y":"N");
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getIntIntnONtoOFF());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getSdIntIntnONtoOFF());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getDutyCyleFromHistogram());
//                psIntegrIntens.print("\t");
//                psIntegrIntens.print(qd.getOnOffThreshold());
////                psIntegrIntens.print(qd.getNumberOfSwitchingEvents());
////                psIntegrIntens.print("\t");
////                psIntegrIntens.print(qd.getNumberOfOnToOffSwitchingEvents());
////                psIntegrIntens.print("\t");
////                psIntegrIntens.print(qd.getNumberOfOffToOnSwitchingEvents());
////                psIntegrIntens.print("\t");
////                psIntegrIntens.print(qd.getSwitchingFrequency());
////                psIntegrIntens.print("\t");
////                psIntegrIntens.print(qd.getOnToOffFrequency());
////                psIntegrIntens.print("\t");
////                psIntegrIntens.print(qd.getOffToOnFrequency());
////                psIntegrIntens.print("\t");
////                psIntegrIntens.print(qd.getDutyCycle());
//                psIntegrIntens.print("\n");
//                IJ.showProgress((double) (dotNo+1)/nDots);
//            }
//            psIntegrIntens.close();
//        } catch(Exception e) {
//            bugFinder.append("Exception in opening and saving file\n");
//        }
//
//    }

    void displayOverlay(ArrayList<MyQuantumDot> dotList) {
        Overlay dotMapOverlay = new Overlay();
        Roi currentRoi = null;
        RoiManager currentRoiManager = new RoiManager(); //define a new Roi manager
        currentRoiManager.reset(); //reset the currentRoiManager
        for(int i=0; i<nDots; i++) {
            currentRoi = new Roi(dotList.get(i).getEnvelope());
            currentRoiManager.addRoi(currentRoi);
            dotMapOverlay.add(currentRoi);
        }

        dotMapOverlay.drawLabels(true); //draw ROIs labels
        imp.setOverlay(dotMapOverlay);
    }

    int[][] readDotMap() {
        //open DotMap file and fill the array int [][] dotMap
        OpenDialog od = new OpenDialog("Select the dotMap file:", "");
        if(od.getFileName()==null) {
            return null;
        }
        File dmFile = new File(od.getPath());
        //find the number of QDs
        int entries = (int) getNentries(dmFile)/2; //two columns
        this.nDots = entries;
        bugFinder.append("Number of entries (dots): "+entries+"\n");

        //create the dotMap array: 2 columns (X and Y coordinates)
        int[][] dotMap = new int[entries][2];
        //exception handling
        try {
            IJ.showStatus("Loading dot map...");
            IJ.showProgress(0.0);
            Reader dmReader = new BufferedReader(new FileReader(dmFile));
            StreamTokenizer dmTok = new StreamTokenizer(dmReader);
            dmTok.eolIsSignificant(false);
            dmTok.parseNumbers();
            int currentEntry=0; //first dot
            boolean isX=true; //fist coordinate is X
            //read the file
            while(dmTok.nextToken()!=dmTok.TT_EOF) {
                if(dmTok.ttype!=dmTok.TT_NUMBER){
                    bugFinder.append("Error in dot map (dot "+(currentEntry+1)+") : token not a number!\n");
                } else {
                    //save the coordinate
                    dotMap[currentEntry][(isX?0:1)]=(int)dmTok.nval;
                    //if it is Y coord. then next time we will have new QD
                    if(!isX) {
                        currentEntry++;
                    }
                    isX = !isX; //to alternate between X and Y coord.
                }
                IJ.showProgress(((double) currentEntry)/entries);
            }
            dmReader.close();
        }catch (IOException e) {
            //meassage
            bugFinder.append("Exception while reading dotMap" + e.getMessage() + "\n");
        }
	// returns an array of int
        return dotMap;
    }

    int[][] readBgDotMap() {
        //open DotMap file and fill the array int bgDotMap[][]
        OpenDialog od = new OpenDialog("Select the background dot map file", "");
        if(od.getFileName()==null) {
            return null;
        }
        File dmFile = new File(od.getPath());
        int entries = (int) getNentries(dmFile)/2; //two columns
        this.nBgDots = entries;
        bugFinder.append("Number of entries (BG dots): "+entries+"\n\n");

        //create bgDotMap: 2 columns (X and Y coordinates)
        int[][] bgDotMap = new int[entries][2];
        //exception handling
        try {
            IJ.showStatus("Loading background dot map...");
            IJ.showProgress(0.0);
            Reader dmReader = new BufferedReader(new FileReader(dmFile));
            StreamTokenizer dmTok = new StreamTokenizer(dmReader);
            dmTok.eolIsSignificant(false);
            dmTok.parseNumbers();
            int currentEntry=0; //first BGdot
            boolean isX=true; //first coordinate (X)
            while(dmTok.nextToken()!=dmTok.TT_EOF) {
                if(dmTok.ttype!=dmTok.TT_NUMBER){
                    bugFinder.append("Error in bg dot map (bg dot "+(currentEntry+1)+") : token not a number!\n");
                } else {
                    bgDotMap[currentEntry][(isX?0:1)]=(int)dmTok.nval;
                    //if it is Y coord. then then next will be a new BGdot
                    if(!isX) {
                        currentEntry++;
                    }
                    isX = !isX;
                }
                IJ.showProgress(((double)currentEntry)/entries);
            }
            dmReader.close();
        }catch (IOException e) {
            //message
            bugFinder.append("Caught IOException while reading bgDotMap: " + e.getMessage() + "\n");
        }
        //returnd an array of int
        return bgDotMap;
    }

    int getNentries(File dmFile){

        bugFinder.append("Getting number of entries from file... \n");
        //fields
        int entries =0;
        //exception handler and calculation of number of entries
        try {
            Reader dmReader = new BufferedReader(new FileReader(dmFile));
            StreamTokenizer dmTok = new StreamTokenizer(dmReader);
            dmTok.eolIsSignificant(false);
            dmTok.parseNumbers();
            //scan all the file
            while(dmTok.nextToken()!=dmTok.TT_EOF) {
                if(dmTok.ttype!=dmTok.TT_NUMBER){
                    bugFinder.append("Error (number of entries in file): token not a number!\n");
                } else {
                    //bugFinder.append("Number: "+(int)dmTok.nval+"\n");
                    entries++;
                }
            }
            dmReader.close();
        }catch (IOException e) {
            //message
            bugFinder.append("Caught IOException while finding entries: " + e.getMessage() + "\n");
       }
        return entries; //each row should have X and Y coordinate
    }

}
