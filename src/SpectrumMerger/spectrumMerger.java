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

        for (Peak p : list.getSortedPeaks(true)) {
            double intensity = p.getIntensity();
            int M2z = (int) Math.round(p.getM2z() * M2Z_ROUNDING_THRESHOLD);

            //System.out.println("INTENSITY: " + intensity + ",  MAXINTENSITY: " + maxIntensity);
            //System.out.println("INTENSITY Normalized: " + (int)(intensity/maxIntensity*INTENSITY_ROUNDING_THRESHOLD));

            double normalizedIntensity = intensity/maxIntensity;
            specArray[M2z] += normalizedIntensity;
        }
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

        for (Complex c : comp){
            //System.out.println(c);
        }

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
