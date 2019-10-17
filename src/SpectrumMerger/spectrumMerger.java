package SpectrumMerger;

import edu.scripps.pms.census.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;

import java.io.IOException;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.*;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

public class spectrumMerger {

    private final static String MS2_FILE_PATH = "/home/ty/timstof/input/data3/N2-pPro-C-01.ms2";
    private final static String MS2_FILE_FORMAT = "ms2";

    private final static int INTENSITY_ROUNDING_THRESHOLD = 1000;
    private final static int M2Z_ROUNDING_THRESHOLD = 1000;
    private final static int SPEC_ARRAY_SIZE = 2097152;

    private double[] specArray;

    private spectrumMerger(Iterator<PeakList> spectraItr) {

        // Create Array with a size of m2z spectrum
        //TODO: Convert to Short Array for better size
        specArray = new double[SPEC_ARRAY_SIZE];

        // combine peaks into one array
        combinePeaks(spectraItr);
    }

    private spectrumMerger() {
        spectrumDbUtil specDB = new spectrumDbUtil();

        if(specDB.connect() == true) {
            specArray = new double[SPEC_ARRAY_SIZE];
            PeakList list = specDB.getSpectrumPeakIntensityByProteinID(1);
            combinePeaks(list);
            specDB.disconnect();
        }
    }

    /*
    Inserts and sums a peaks intensity into the combined spectra array by the index: mass*1000
    */
    private void combinePeaks(Iterator<PeakList> spectraItr) {

        while (spectraItr.hasNext()) {
            PeakList list = spectraItr.next();

            // used to normalize the peakList
            double maxIntensity = list.getMaxIntensity();

            for (Peak p : list.getSortedPeaks(true)) {
                double intensity = p.getIntensity();
                int M2z = (int) Math.round(p.getM2z() * M2Z_ROUNDING_THRESHOLD);

                //System.out.println("INTENSITY: " + intensity + ",  MAXINTENSITY: " + maxIntensity);
                //System.out.println("INTENSITY Normalized: " + (int)(intensity/maxIntensity*INTENSITY_ROUNDING_THRESHOLD));

                //                specArray[M2z] += (int)(intensity/maxIntensity*INTENSITY_ROUNDING_THRESHOLD);

                double normalizedIntensity = intensity/maxIntensity;
                specArray[M2z] += normalizedIntensity;
            }
        }
    }

    /*
    Inserts and sums a peaks intensity into the combined spectra array by the index: mass*1000
    */
    private void combinePeaks(PeakList list) {

        // used to normalize the peakList
        double maxIntensity = list.getMaxIntensity();

        int count = 0;

        for (Peak p : list.getSortedPeaks(true)) {
            count++;
            double intensity = p.getIntensity();
            int M2z = (int) Math.round(p.getM2z() * M2Z_ROUNDING_THRESHOLD);

            //System.out.println("INTENSITY: " + intensity + ",  MAXINTENSITY: " + maxIntensity);
            //System.out.println("INTENSITY Normalized: " + (int)(intensity/maxIntensity*INTENSITY_ROUNDING_THRESHOLD));

            double normalizedIntensity = intensity/maxIntensity;
            specArray[M2z] += normalizedIntensity;
        }

        System.out.println(count);
    }

    private ArrayList<Peak> normalize(ArrayList<Peak> peakList){

        double maxIntensity = 0;

        for(int i = 0; i < peakList.size(); i++){
            if(peakList.get(i).getIntensity() > maxIntensity){
                maxIntensity = peakList.get(i).getIntensity();
            }
        }

        return null;
    }

    /*
    Finds Minimum Mass value for spectrumArray
    */
    private double getMinM2z() {
        return convertToMass(getMinM2zIndex());

    }

    /*
    Finds Minimum Mass index for spectrumArray
    */
    private int getMinM2zIndex() {
        for(int i = 0; i < specArray.length; i++){
            if(specArray[i] != 0) return i;
        }
        return -1;
    }

    /*
    Finds Maximum Mass value for spectrumArray
    */
    private double getMaxM2z() {
        return convertToMass(getMaxM2zIndex());
    }

    /*
    Finds Maximum Mass index for spectrumArray
    */
    private int getMaxM2zIndex() {
        for(int i = specArray.length -1; i > 0; i--){
            if(specArray[i] != 0) return i;
        }
        return -1;
    }

    private Complex[] ftrans(){

        for (int i = 0; i<specArray.length; i++){
            //System.out.println(specArray[i]);
        }
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] comp = fft.transform(specArray, TransformType.FORWARD);

        return comp;
    }

    /*
    Converts specArray index to its original m2z data
    */
    private double convertToMass(int i){
        return ((double)i)/M2Z_ROUNDING_THRESHOLD;
    }

    private double getIntesity(){
    return 1.0;
    }

    private short setIntesity(){
        return (short)(1);
    }

    private int getNumberPeaks(){
        int counter = 0;
        for(double intensity : specArray){
            if(intensity != 0.0) {
                //System.out.println(intensity);
                counter++;
            }
        }
        return counter;
    }

    private int getNumberPeaks(Complex[] comp){
        int counter = 0;
        for(Complex c : comp){
            if(c.getImaginary() != 0.0 || c.getReal() != 0.0) counter++;
        }
        return counter;
    }

    public void printCompData(Complex[] comp){

        double meanReal = 0;
        double meanImj = 0;
        double meanComb = 0;


        double minReal = Double.MAX_VALUE;
        double maxReal = Double.MIN_VALUE;
        double minImj = Double.MAX_VALUE;
        double maxImj = Double.MIN_VALUE;

        double minComb = Double.MAX_VALUE;
        double maxComb = Double.MIN_VALUE;

        double sdReal = 0;
        double sdImj = 0;
        double sdComb = 0;

        int pointsAbove20Real = 0;
        int pointsAbove20Imj = 0;
        int pointsAbove20Comb = 0;

        int pointsAbove10Real = 0;
        int pointsAbove10Imj = 0;
        int pointsAbove10Comb = 0;

        int pointsAbove30Comb = 0;
        int pointsAbove40Comb = 0;




        double[] peaks = new double[20];


        int counter = 1;
        for(Complex c : comp){

            double dist = Math.pow(c.getReal(),2) + Math.pow(c.getImaginary(), 2);
            dist = Math.sqrt(dist);

            meanReal = (meanReal + c.getReal())/counter;
            meanImj = (meanImj + c.getImaginary())/counter;
            meanComb = (meanComb + dist)/counter;

            if(minReal >= c.getReal()) minReal = c.getReal();
            if(maxReal <= c.getReal()) maxReal = c.getReal();

            if(minImj >= c.getImaginary()) minImj = c.getImaginary();
            if(maxImj <= c.getImaginary()) maxImj = c.getImaginary();


            if(minComb >= dist) minComb = dist;
            if(maxComb <= dist) maxComb = dist;

            if(c.getReal() > 20) {
                pointsAbove20Real++;
            }
            if(c.getImaginary() > 20) {
                pointsAbove20Imj++;
            }

            if(c.getImaginary() > 20) {
                pointsAbove20Imj++;
            }
            if(dist > 20) pointsAbove20Comb++;

            if(c.getReal() > 10) pointsAbove10Real++;
            if(c.getImaginary() > 10) pointsAbove10Imj++;
            if(dist > 10) pointsAbove10Comb++;

            if(dist > 30) pointsAbove30Comb++;
            if(dist > 67) pointsAbove40Comb++;



            counter++;
        }

        counter = 1;
        for(Complex c : comp){
            double dist = Math.pow(c.getReal(),2) + Math.pow(c.getImaginary(), 2);

            sdReal += Math.pow((c.getReal() - meanReal),2) / comp.length;
            sdImj += Math.pow((c.getImaginary() - meanImj),2) / comp.length;
            sdComb += Math.pow((dist - meanComb),2) / comp.length;

            counter++;
        }


        System.out.println("--------mean--------");
        //System.out.println("meanReal: " + meanReal);
        //System.out.println("meanImj: " + meanImj);
        System.out.println("meanComb: " + meanComb);

        System.out.println("--------mix/max--------");
        //System.out.println("minReal: " + minReal);
        //System.out.println("maxReal: " + maxReal);
        //System.out.println("minImj: " + minImj);
        //System.out.println("maxImj: " + maxImj);
        System.out.println("minComb: " + minComb);
        System.out.println("maxComb: " + maxComb);

        System.out.println("--------sd--------");
        //System.out.println("sdReal: " + sdReal);
        //System.out.println("sdImj: " + sdImj);
        System.out.println("sdComb: " + sdComb);

        System.out.println("--------points above 20--------");
        //System.out.println("pointsAbove20Real: " + pointsAbove20Real);
        //System.out.println("pointsAbove20Imj: " + pointsAbove20Imj);
        System.out.println("pointsAbove20Comb: " + pointsAbove20Comb);

        System.out.println("--------points above 10--------");
        //System.out.println("pointsAbove10Real: " + pointsAbove10Real);
        //System.out.println("pointsAbove10Imj: " + pointsAbove10Imj);
        System.out.println("pointsAbove10Comb: " + pointsAbove10Comb);

        System.out.println("pointsAbove30Comb: " + pointsAbove30Comb);
        System.out.println("pointsAbove40Comb: " + pointsAbove40Comb);


    }

    public double[] getSpecArr(){
        return specArray;
    }


    public static void main(String[] args) throws IOException, SQLException {

        spectrumMerger sm = new spectrumMerger();

        //sm.printStats();
        //sm.printSpectra();
        Complex[] comp = sm.ftrans();

        System.out.println(comp.length);
        System.out.println(sm.specArray.length);

        System.out.println(comp.length);
        System.out.println(sm.specArray.length);
        System.out.println(sm.getNumberPeaks(comp));
        System.out.println(sm.getNumberPeaks());

        //sm.printCompData(comp);
        fileManager fm = new fileManager();
        fm.writeFFTtoCSV(comp);
        fm.writeCombinedSpectrumtoCSV(sm.getSpecArr());
        //for (Complex c : comp){
            //System.out.println(c);
        //}

        /*
        SpectrumReader sr = new SpectrumReader(MS2_FILE_PATH, MS2_FILE_FORMAT);
        Iterator<PeakList> spectraItr = sr.getSpectra();

        spectrumMerger sm = new spectrumMerger(spectraItr);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);




        //fft.transform(specArray,DftNormalization.STANDARD);
        //sm.ftrans();
        //sm.printStats();
        //sm.printSpectra();

        sr.closeDataFile();
        System.out.println("Finished");
        */
    }



    ////////////////////////////////////PRINTS/////////////////////////////////////////////

    /*
    For Tests
    */
    private void printStats(){
        System.out.println("Max Mass Value: " + getMaxM2z());
        System.out.println("Max Mass Index: " + getMaxM2zIndex());
        System.out.println("Min Mass Value: " + getMinM2z());
        System.out.println("Min Mass Index: " + getMinM2zIndex());
    }

    /*
    Splits the combined spectra into 200 buckets of equal size mass and prints the average of the buckets intensities
    */
    private void printSpectra(){

        int spectra200BucketThresh = 10000;

        // summedRepresentation = condensered version of specArray, with 200 buckets
        double[] summedRepresentation = new double[SPEC_ARRAY_SIZE/spectra200BucketThresh];

        // Sum every 10000 elements from specArray into one bucket of summedRepresentation
        for(int i = 0; i<summedRepresentation.length; i++){
            double sum = 0;
            for(int j = i*spectra200BucketThresh; j < (i+1)*spectra200BucketThresh; j++){
                sum += specArray[j];
            }
            summedRepresentation[i] = sum/spectra200BucketThresh;
        }

        System.out.println("summedRepresentation length: " + summedRepresentation.length);

        int counter = 0;
        for(double sum : summedRepresentation){
            StringBuilder s = new StringBuilder();
            for(int i = 0; i<sum; i++){
                s.append(".");
            }
            System.out.println("Mass: [ " + counter*spectra200BucketThresh/M2Z_ROUNDING_THRESHOLD + " - "
                    + (counter+1)*spectra200BucketThresh/M2Z_ROUNDING_THRESHOLD + " ] : " + sum + "  | " + s );
            counter++;
        }
    }

}
