package SpectrumMerger;

import edu.scripps.pms.util.spectrum.PeakList;
import javafx.util.Pair;

import java.io.IOException;

class SpectrumMergerHandler {


    public static void main(String[] args) throws IOException {
        generateAllPeptideSpectrum(0.005);
    }

    public static void generateAllPeptideSpectrum(double threshold) throws IOException {

        int numOfPeptideIds = 10;
        for(int peptideID = 0; peptideID < numOfPeptideIds; peptideID++){

            fileManager fm = new fileManager();
            spectrumDbUtil specDB = new spectrumDbUtil();
            spectrumMerger sm =  new spectrumMerger(peptideID, threshold);

            PeakList[] individualPeptidePeakLists;
            if(specDB.connect() == true) {
                individualPeptidePeakLists = specDB.getSpectrumPeakIntensityByProteinID2(peptideID);
                specDB.disconnect();
            }else{
                System.out.println("Error Opening connection to Database");
                return;
            }
            sm.combinePeaks(individualPeptidePeakLists);
            sm.mergePeaksByThreshold();

            fm.writeCombinedSpectrumtoCSV(sm.getSpecArr(), peptideID);

            System.out.println("Done");
        }

    }

    public static void generateSinglePeptideSpectrum(int peptideID, double threshold) throws IOException {

        fileManager fm = new fileManager();
        spectrumDbUtil specDB = new spectrumDbUtil();
        spectrumMerger sm =  new spectrumMerger(peptideID, threshold);

        PeakList[] individualPeptidePeakLists;
        if(specDB.connect() == true) {
            individualPeptidePeakLists = specDB.getSpectrumPeakIntensityByProteinID2(peptideID);
            specDB.disconnect();
        }else{
            System.out.println("Error Opening connection to Database");
            return;
        }
        sm.combinePeaks(individualPeptidePeakLists);
        sm.mergePeaksByThreshold();

        fm.writeCombinedSpectrumtoCSV(sm.getSpecArr(), peptideID);

        System.out.println("Done");

    }

    public static void generateSingleIndividualCSV(int peptideId, double threshold) throws IOException {
        fileManager fm = new fileManager();
        spectrumDbUtil specDB = new spectrumDbUtil();
        spectrumMerger sm =  new spectrumMerger(peptideId, threshold);

        PeakList[] individualPeptidePeakLists;

        // get array containing peak lists of each peptide sample
        if(specDB.connect() == true) {
            individualPeptidePeakLists = specDB.getSpectrumPeakIntensityByProteinID2(1);
            System.out.println("Number of peaks: " + individualPeptidePeakLists.length);

            specDB.disconnect();
        }else{
            return;
        }

        //Generate sample individual graphs
        fm.writeIndividualSpectrumtoCSV(individualPeptidePeakLists, 10);
    }

}
