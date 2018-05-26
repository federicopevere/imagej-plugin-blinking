/*
* class MyQuantumDot: it represent a quantum dot in a photoluminescence blinking trace
*
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

//import packages
import ij.*;
import ij.measure.*;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.stat.StatUtils;

//class MyQuantumDot: it represent a quantum dot in a PL blinking trace
public class MyQuantumDot {

    //attributes of the dot
    int centerX = 0; //center of the QD: x coordinate
    int centerY = 0;//center of the QD: y coordinate
    Rectangle envelope = null; //rectangle selection
    double[] centerOfMassX = null; //center of mass through sequence
    double[] centerOfMassY = null; //center of mass through sequence
    int diameter = 0; //QD rectangle edge length in pixels
    int bgDistance=0; //distance of bg pixels from edge of the envelope
    int slices = 0; //number of slices of image stack
    double[] integrIntens = null; // sum of intensities of pixels included in QD selection area
    double[] bgIntegrIntens = null; // average background integrated intensity per frame for a quantum dot
    double[] bgCorrIntegrIntens = null;
    double[] maxIntens = null; //maximum intensity of pixels included in dot in a frame
    double[] minIntens = null; //minimum intensity of pixels included in dot in a frame

    //statistical attributes, useful for plotting and making histograms
    double maxIntegrIntens = Double.MIN_VALUE;
    double minIntegrIntens = Double.MAX_VALUE;
    double meanIntegrIntens = 0.0;
    double stddevIntegrIntens = 0.0;

    double maxBgIntegrIntens = Double.MIN_VALUE;
    double minBgIntegrIntens = Double.MAX_VALUE;
    double meanBgIntegrIntens = 0.0;
    double stddevBgIntegrIntens = 0.0;

    double meanBgCorrIntegrIntens=0.0;
    double maxBgCorrIntegrIntens=Double.MIN_VALUE;
    double minBgCorrIntegrIntens=Double.MAX_VALUE;
    double stddevBgCorrIntegrIntens=0.0;

    //attributes which can be evaluated from the integr. intens. histogram
    double meanIntegrIntensON = 0;
    double areaON = 0;
    double stddevIntegrIntensON = 0;
    double meanIntegrIntensOFF = 0;
    double areaOFF = 0;
    double stddevIntegrIntensOFF = 0;
    double dutyCycleHist = 0;
    double lowerThr=0;
    double upperThr=0;

    //attributes which can be evaluated from the blinking trace and
    //the thresholds
    int switchingEvents = 0;
    ArrayList<Integer> onTimes = new ArrayList<Integer>();
    ArrayList<Integer> offTimes = new ArrayList<Integer>();
    double meanOnTime = 0;
    double meanOffTime = 0;
    ArrayList<Integer> onTimesFreqCounts = new ArrayList<Integer>();
    ArrayList<Integer> offTimesFreqCounts = new ArrayList<Integer>();
    ArrayList<Integer> onTimesProbDens = new ArrayList<Integer>();
    ArrayList<Integer> offTimesProbDens = new ArrayList<Integer>();
    double onTimeConstant = 0;
    double offTimeConstant = 0;
    double dutyCycleStat = 0;

    //class constructor: it calculates  the envelope (rectangle) and it extracts
    //the integrated intensity from the blinking trace
    public MyQuantumDot(int centerX, int centerY, int diameter, int bgDistance, ImagePlus imp){
        //assign attributes of the QD
        this.centerX = centerX; //assign X coord
        this.centerY = centerY; //assign Y coord
        this.diameter = diameter; //assign diameter
        this.slices = imp.getNSlices(); //assign n. of slices or frames
        this.bgDistance = bgDistance; //assign bg distance
        //determine the envelope
        determineEnvelope(centerX, centerY, diameter, imp.getProcessor());
        //extract all the integrated intensities (signal+bg, bg, signal)
        extractAllIntegrIntens(imp);
        //calculate statistical information of all the integrated intensities
        //(signal+bg, bg, signal)
        calculateStatAllIntegrIntens();
    }

    //method that creates the envelope (rectangle object) for the dot
    //NB: in Rectangle class X and Y represent the upper left-most corner
    private void determineEnvelope(int centerX, int centerY, int diameter, ImageProcessor ip) {
        //even number of pixels: the coordinates (centerX, centerY) indicate
        //one of the four pixels present at the center of the QD. To find the
        //correct pixel, we calculated the integr. intensity in those four cases
        if(diameter%2==0){
            //integrate intensity for a rectangle having a bigger diameter (1 px
            //difference), i.e. of size (diameter+1)*(diameter+1)
            int tempDiam = diameter+1;
            int[] tempBlock = new int[tempDiam];
            for(int y=0; y<tempDiam; y++) {
                for(int x=0; x<tempDiam; x++){
                    tempBlock[y*tempDiam+x] = ip.getPixel(x+centerX-(int)((tempDiam-1)/2), y+centerY-(int)((tempDiam-1)/2));
                }
            }
            //caclulate integr. intens. of each rectangle of size
            //diameter*diameter
            int[] sum = {0,0,0,0};
            for(int y = 0; y<diameter; y++) {
                for(int x=0; x<diameter; x++) {
                    sum[0] += tempBlock[y*tempDiam+x]; // top left rectangle
                    sum[1] += tempBlock[y*tempDiam+x+1]; // top right rect.
                    sum[2] += tempBlock[(y+1)*tempDiam+x]; // bottom left rect.
                    sum[3] += tempBlock[(y+1)*tempDiam+x+1]; // bottom right rect.
                }
            }
            int maxIndex = 0;
            int maxVal = sum[0];
            for(int i = 1; i<4; i++) {
                if(sum[i]>maxVal) {
                    maxVal = sum[i];
                    maxIndex = i;
                }
            }
            this.envelope = new Rectangle(centerX-(int)diameter/2-((maxIndex==0 || maxIndex==2)?1:0),centerY-(int)diameter/2-((maxIndex==0 || maxIndex==1)?1:0), diameter, diameter);
        } else {
            this.envelope = new Rectangle(centerX-(int)(diameter-1)/2, centerY-(int)(diameter-1)/2, diameter, diameter);
        }
    }

    //FILL: integrIntens[], maxIntens[], minIntens[], maxIntegrIntens,
    //      minIntegrIntens, centerOfMassX[], centerOfMassY[],
    //      bgIntegrIntens[], bgCorrIntens[]
    private void extractAllIntegrIntens(ImagePlus imp) {
        //go through slides and use envelope to extract pixels
        //sum pixels and save sum
        this.integrIntens = new double[this.slices];
        this.maxIntens = new double[this.slices];
        this.minIntens = new double[this.slices];
        this.centerOfMassX = new double[this.slices];
        this.centerOfMassY = new double[this.slices];
        this.bgIntegrIntens = new double[this.slices];
        this.bgCorrIntegrIntens = new double[this.slices];
        ImageStack stack = imp.getStack();
        ImageProcessor ip = null;
        int eX = (int) this.envelope.getX(); //envelope X (left-most X coord.)
        int eY = (int) this.envelope.getY(); //envelope Y (top Y coord.)
        int width = (int) this.envelope.getWidth(); //envelope width
        int height = (int) this.envelope.getHeight(); //envelope height
        int imageWidth = imp.getWidth();
        int imageHeight = imp.getHeight();
        int[] pixels = new int[width*height];
        int bgDist=this.bgDistance;

        //NB the slices are numbered from 1 so the loop index must start from 1
        for(int currentFrame=1; currentFrame<=this.slices; currentFrame++) {

            ip = stack.getProcessor(currentFrame); //get current frame

            //initialization
            int maxIntPx = Integer.MIN_VALUE;
            int minIntPx = Integer.MAX_VALUE;
            int sumCenterOfMassX = 0;
            int sumCenterOfMassY = 0;
            int bgPixels = 0;
            int tempBGsum = 0;
            int sum = 0;
            int p = 0;

            //loop inside the envelope: we start from the top left-most corner
            //of the envelope (eX,eY)
            for(int y=eY; y<eY+height; y++) {
                for(int x=eX; x<eX+width; x++) {

                    pixels[p] = ip.getPixel(x,y);
                    //calculate center of mass
                    sumCenterOfMassX += x*pixels[p];
                    sumCenterOfMassY += y*pixels[p];
                    //calculate background:
                    //out of top left-most corner (inside image)
                    if(x==eX && y==eY && eX>(bgDist-1) && eY>(bgDist-1)) {
                        tempBGsum += ip.getPixel(x-bgDist,y-bgDist);
                        bgPixels++;
                    }
                    //out of bottom left-most corner (inside image)
                    if(x==eX && y==eY+height-1 && eX>(bgDist-1) && eY<(imageHeight-bgDist)) {
                        tempBGsum += ip.getPixel(x-bgDist,y+bgDist);
                        bgPixels++;
                    }
                    //out of top right-most corner (inside image)
                    if(x==eX+width-1 && y==eY && eX<(imageWidth-bgDist) && eY>(bgDist-1)) {
                        tempBGsum += ip.getPixel(x+bgDist,y-bgDist);
                        bgPixels++;
                    }
                    //out of bottom right-most corner (inside image)
                    if(x==eX+width-1 && y==eY+height-1 && eX<(imageWidth-bgDist) && eY<(imageHeight-bgDist)) {
                        tempBGsum += ip.getPixel(x+bgDist,y+bgDist);
                        bgPixels++;
                    }
                    //left edge (inside image)
                    if(x==eX && eX>0) {
                        tempBGsum += ip.getPixel(x-bgDist,y);
                        bgPixels++;
                    }
                    //top edge (inside image)
                    if(y==eY && eY>0) {
                        tempBGsum += ip.getPixel(x,y-bgDist);
                        bgPixels++;
                    }
                    //right edge (inside image)
                    if(x==eX+width-1 && x<imageWidth-1) {
                        tempBGsum += ip.getPixel(x+bgDist,y);
                        bgPixels++;
                    }
                    //bottom edge (inside image)
                    if(y==eY+height-1 && eY<imageHeight-1) {
                        tempBGsum += ip.getPixel(x,y+bgDist);
                        bgPixels++;
                    }
                    //max intensity pixel
                    if(pixels[p]>maxIntPx){
                        maxIntPx = pixels[p];
                    }
                    //min intensity pixel
                    if(pixels[p]<minIntPx) {
                        minIntPx = pixels[p];
                    }

                    sum += pixels[p]; //integrated intensity
                    p++; //next pixel
                }
            }
            //save the values in the related fields of the MyQuantumDot class
            this.centerOfMassX[currentFrame-1] = (double) sumCenterOfMassX/sum;
            this.centerOfMassY[currentFrame-1] = (double) sumCenterOfMassY/sum;
            this.maxIntens[currentFrame-1] = maxIntPx;
            this.minIntens[currentFrame-1] = minIntPx;
            this.integrIntens[currentFrame-1] = sum;
            this.bgIntegrIntens[currentFrame-1] = tempBGsum/bgPixels*diameter*diameter;
            this.bgCorrIntegrIntens[currentFrame-1] = this.integrIntens[currentFrame-1]-this.bgIntegrIntens[currentFrame-1];
        }
    } //end of class extractIntegrIntens

    private void calculateStatAllIntegrIntens(){
        //statistical data, useful for plotting and making histograms
        this.maxIntegrIntens = StatUtils.max(this.integrIntens);
        this.minIntegrIntens = StatUtils.min(this.integrIntens);
        this.meanIntegrIntens = StatUtils.mean(this.integrIntens);
        this.stddevIntegrIntens = Math.sqrt(StatUtils.variance(this.integrIntens));

        this.maxBgIntegrIntens = StatUtils.max(this.bgIntegrIntens);
        this.minBgIntegrIntens = StatUtils.min(this.bgIntegrIntens);
        this.meanBgIntegrIntens = StatUtils.mean(this.bgIntegrIntens);
        this.stddevBgIntegrIntens = Math.sqrt(StatUtils.variance(this.bgIntegrIntens));

        this.maxBgCorrIntegrIntens = StatUtils.max(this.bgCorrIntegrIntens);
        this.minBgCorrIntegrIntens = StatUtils.min(this.bgCorrIntegrIntens);
        this.meanBgCorrIntegrIntens = StatUtils.mean(this.bgCorrIntegrIntens);
        this.stddevBgCorrIntegrIntens = Math.sqrt(StatUtils.variance(this.bgCorrIntegrIntens));
    }

    public int getX() {
        return this.centerX;
    }

    public int getY() {
        return this.centerY;
    }

    public double getIntegrIntens(int frame) {
        return this.integrIntens[frame];
    }

    public double[] getIntegrIntens() {
        return this.integrIntens;
    }

    public double getBgIntegrIntens(int frame) {
        return this.bgIntegrIntens[frame];
    }

    public double[] getBgIntegrIntens() {
        return this.bgIntegrIntens;
    }

    public double getBgCorrIntegrIntens(int frame) {
        return this.bgCorrIntegrIntens[frame];
    }

    public double[] getBgCorrIntegrIntens() {
        return this.bgCorrIntegrIntens;
    }

    //statistical quantities

    public double getMeanIntegrIntens(){
        return this.meanIntegrIntens;
    }

    public double getStddevIntegrIntens(){
        return this.stddevIntegrIntens;
    }

    public double getMaxIntegrIntens(){
        return this.maxIntegrIntens;
    }

    public double getMinIntegrIntens(){
        return this.minIntegrIntens;
    }

    public double getMeanBgIntegrIntens(){
        return this.meanBgIntegrIntens;
    }

    public double getStddevBgIntegrIntens(){
        return this.stddevBgIntegrIntens;
    }

    public double getMaxBgIntegrIntens(){
        return this.maxBgIntegrIntens;
    }

    public double getMinBgIntegrIntens(){
        return this.minBgIntegrIntens;
    }

    public double getMeanBgCorrIntegrIntens(){
        return this.meanBgCorrIntegrIntens;
    }

    public double getStddevBgCorrIntegrIntens(){
        return this.stddevBgCorrIntegrIntens;
    }

    public double getMaxBgCorrIntegrIntens(){
        return this.maxBgCorrIntegrIntens;
    }

    public double getMinBgCorrIntegrIntens(){
        return this.minBgCorrIntegrIntens;
    }

    public Rectangle getEnvelope() {
        return this.envelope;
    }

    public double getMeanIntegrIntensON(){
        return this.meanIntegrIntensON;
    }

    public double getStddevIntegrIntensON(){
        return this.stddevIntegrIntensON;
    }

    public double getMeanIntegrIntensOFF(){
        return this.meanIntegrIntensOFF;
    }

    public double getStddevIntegrIntensOFF(){
        return this.stddevIntegrIntensOFF;
    }

    public int getSwitchingEvents(){
        return this.switchingEvents;
    }

    public double getDutyCycleHist(){
        return this.dutyCycleHist;
    }

    public void setBgIntegrIntens(int frame, int value) {
        this.bgIntegrIntens[frame] = value;
    }

    //incomplete and unused at the moment
    public void recalculateIntegrIntens(ImageProcessor ip, int frame, int eX, int eY) {
        if(this.diameter%2==0) {
            //determine if local envelope has to be readjusted
        }
        int width = (int) this.envelope.getWidth();
        int height = (int) this.envelope.getHeight();
        int sum = 0;
        for(int y=eY; y<eY+height; y++) {
            for(int x=eX; x<eX+width; x++) {
                sum += ip.getPixel(x,y);
            }
        }
        this.integrIntens[frame-1] = sum;
    }

} //end of class MyQuantumDot
