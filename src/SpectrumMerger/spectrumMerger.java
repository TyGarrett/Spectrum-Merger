package SpectrumMerger;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;

import java.io.IOException;
import java.util.Iterator;

import javafx.util.Pair;
import org.apache.commons.math3.complex.Complex;

import java.sql.SQLException;
import java.util.LinkedList;

public class spectrumMerger {

    private final static int M2Z_ROUNDING_THRESHOLD = 1000;
    private final static int SPEC_ARRAY_SIZE = 2097152;
    private final static int PPM_DIVISOR = 1000000;

    // Index is mass*M2Z_ROUNDING_THRESHOLD.
    // Value is the intensity
    private float[] specArray;

    // Index represents the intensity of peaks in specArray[].
    // Value at each index represents the mass in specArray[].
    private int[] sortedSpecArrayByIndex;
    private int numberOfNonZeroElements;
    private int peptideID;

    private int massKeyRunningAverage = 0;
    private int massKeyNumberOfItems = 0;


    /*
    Constructor for spectrumMerger. Creates empty specArray.
    */
    spectrumMerger(int peptideId) {
        specArray = new float[SPEC_ARRAY_SIZE];
        this.peptideID = peptideId;
    }

    /*
    Inserts and sums a peaks intensity into the combined spectra array by the index: mass*1000
    */
    private void addPeakListToSpecArray(PeakList list) {

        // used to normalize the peakList
        double maxIntensity = list.getMaxIntensity();

        for (Peak p : list.getSortedPeaks(true)) {
            double intensity = p.getIntensity();
            int M2z = (int) Math.round(p.getM2z() * M2Z_ROUNDING_THRESHOLD);
            float normalizedIntensity = (float)(intensity/maxIntensity);
            specArray[M2z] += normalizedIntensity;
        }
    }

    /*
       adds a a list of pairs into the specArray
     */
    void addPairToSpecArray(LinkedList<Pair<Float, Float>> list) {
        for (Pair p : list) {
            float intensity = (float) p.getValue();
            int M2z = (int) Math.round(((float) p.getKey()) * M2Z_ROUNDING_THRESHOLD);
            specArray[M2z] += intensity;
        }
    }

    /*
    adds each peak in list to the specArray
    */
    void addPeaksToSpecArray(PeakList[] list){
        for(PeakList pl : list){
            addPeakListToSpecArray(pl);
        }
    }

    /*
    creates an array which contains the sorted index's of specArray by intensity
    */
    private void generateSortedSpecArrayByIndex(){
        sortedSpecArrayByIndex = new int[getNumberPeaks()];

        //fill sorted Array with all index's of spec array
        int sortedIdx = 0;
        for(int i = 0; i < specArray.length; i++){
            if(specArray[i] != 0){
                sortedSpecArrayByIndex[sortedIdx] = i;
                sortedIdx++;
            }
        }
        quickSort(sortedSpecArrayByIndex, 0, sortedSpecArrayByIndex.length-1);
    }

    /**
     * sorts arr[] based on the intensity of specArray[]
     */
    private void quickSort(int[] arr, int low, int high)
    {
        if (low < high)
        {
            /* pi is partitioning index, arr[pi] is
              now at right place */
            int pi = partition(arr, low, high);

            // Recursively sort elements before
            // partition and after partition
            quickSort(arr, low, pi-1);
            quickSort(arr, pi+1, high);
        }
    }

    private int partition(int[] arr, int low, int high)
    {
        int pivot = arr[high];
        int i = (low-1); // index of smaller element
        for (int j=low; j<high; j++){
            // If current element is smaller than or
            // equal to pivot
            if (specArray[arr[j]] <= specArray[pivot]){
                i++;
                // swap arr[i] and arr[j]
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }

        // swap arr[i+1] and arr[high] (or pivot)
        int temp = arr[i+1];
        arr[i+1] = arr[high];
        arr[high] = temp;

        return i+1;
    }

    /**
     * merge peaks together: start with largest and increments to smallest.
     */
    void mergePeaksByThreshold(int ppm){

        generateSortedSpecArrayByIndex();

        for (int i = sortedSpecArrayByIndex.length-1; i >= 0; i--){

            // Find the next Biggest peak. skip already merged peaks
            if(specArray[sortedSpecArrayByIndex[i]] == 0){
                continue;
            }

            // Index of spec array with largest peak
            int peakIdx = sortedSpecArrayByIndex[i];
            int threshold = getThreshold(ppm, peakIdx);
            this.numberOfNonZeroElements++;

            //Add Peaks to teh right of peak that are within the threshold
            for(int j = peakIdx + 1; j <= (peakIdx + threshold) && j < specArray.length; j++){
                specArray[peakIdx] += specArray[j];
                specArray[j] = 0;
            }

            //Add Peaks to the left of threshold that are within the threshold
            for(int j = peakIdx - 1; j >= (peakIdx - threshold) && j >= 0; j--){
                specArray[peakIdx] += specArray[j];
                specArray[j] = 0;
            }
        }
    }

    /**
     * converts spectrumArray into a mz and intensity array
     */
    Pair<float[], float[]> getMassIntensityArrays(){
        float[] peakMzArray = new float[numberOfNonZeroElements];
        float[] peakIntensityArray = new float[numberOfNonZeroElements];

        int mzIndex = 0;
        for(int i = 0; i < specArray.length; i++) {
            if (specArray[i] != 0) {
                peakMzArray[mzIndex] = ((float) i) / 1000;
                peakIntensityArray[mzIndex] = specArray[i];
                mzIndex++;
            }
        }

        return new Pair<float[], float[]>(peakMzArray, peakIntensityArray);
    }

    /**
     * calculates the threshold for which to merge peaks associated PPM
     */
    int getThreshold(int ppm, int peakIdx){
        return ((ppm*peakIdx)/(PPM_DIVISOR));
    }

    int getMassKey(){
        return this.massKeyRunningAverage;
    }

    void setMassKey(int val){
        massKeyNumberOfItems += 1;
        massKeyRunningAverage = (massKeyRunningAverage*(massKeyNumberOfItems-1) + val) / massKeyNumberOfItems;
    }


    public int getNumberOfNonZeroElements(){
        return this.numberOfNonZeroElements;
    }

    /*
    Converts specArray index to its original m2z data
    */
    private double convertToMass(int i){
        return ((double)i)/M2Z_ROUNDING_THRESHOLD;
    }

    private int getNumberPeaks(){
        int counter = 0;
        for(float intensity : specArray){
            if(intensity != 0.0) {
                //System.out.println(intensity);
                counter++;
            }
        }
        return counter;
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


    float[] getSpecArr(){
        return specArray;
    }

    public int getPeptideID(){
        return peptideID;
    }

    //////////////////////////PRINTS/////////////////////////////////////////////

    /*
    For Tests
    */
    private void printStats(){
        System.out.println("Max Mass Value: " + getMaxM2z());
        System.out.println("Max Mass Index: " + getMaxM2zIndex());
        System.out.println("Min Mass Value: " + getMinM2z());
        System.out.println("Min Mass Index: " + getMinM2zIndex());
    }

      /////////////////////////OLD METHODS POSSIBLY DELETE LATER////////////////////////////


}
