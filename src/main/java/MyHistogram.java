/*************************************************************************
 *  Compilation:  javac MyHistogram.java
 *
 *  This data type supports simple client code to create dynamic
 *  histograms of the frequency of occurrence of values in [0, N).
 *  The frequencies are kept in an instance-variable array, and
 *  an instance variable max tracks the maximum frequency (for scaling).
 *
 *  % java MyHistogram 50 1000000 
 *
 *************************************************************************/

public class MyHistogram {
    private int [] freq;   // freq[i] = # occurences of value i
    private double minBinStart; // min bin start
    private double maxBinEnd; // max bin end
    private int freqMax;  // max frequency of any value

    // Create a new histogram. 
    public MyHistogram(int Nbins, double minBinStart, double maxBinEnd) {
        this.freq = new int[Nbins];
        this.minBinStart = minBinStart;
        this.maxBinEnd = maxBinEnd;
    }

    // Add one occurrence of the value i. 
    public void addDataPoint(int i) {
        freq[i]++; 
        if (freq[i] > freqMax) freqMax = freq[i]; 
    }
    
    public double getMinBinStart(){
       return this.minBinStart;
    }
    
    public double getMaxBinEnd(){
       return this.maxBinEnd;
    }
    
    public int[] getFreq(){
       return this.freq;
    }
    
    public int getFreq(int i){
        return this.freq[i];
    }
    
    public double getBinCenter(int i){
        double histogramRange = this.maxBinEnd-this.minBinStart;
        int histogramNbins = this.freq.length;
        double histogramBinSize = histogramRange/histogramNbins;  
        return (this.minBinStart+histogramBinSize/2+i*histogramBinSize);
    }
    
    public double[] getBinCenter(){
        double[] binCenter=null;
        for(int i=0; i<this.freq.length; i++){
            binCenter[i]=this.getBinCenter(i);
        }
        return binCenter;
    }
    
    public double getRange(){
        return this.maxBinEnd-this.minBinStart;
    }
    
    public int getNbins(){
        return this.freq.length;
    }
    
    public double getBinSize(){
        return this.getRange()/((double) this.getNbins());
    }
    
    public int getFreqMax(){
       return this.freqMax;
    }

} 
