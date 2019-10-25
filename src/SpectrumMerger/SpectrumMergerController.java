package SpectrumMerger;

import edu.scripps.pms.util.spectrum.PeakList;
import org.jboss.util.Null;

import java.io.File;
import java.io.IOException;

class SpectrumMergerHandler {
    spectrumDbUtil copyDBMT = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");

    public SpectrumMergerHandler(){

    }

    public static void main(String[] args) throws IOException {
        String databaseFileName = "testLibDuplicateSpectra.db";
        double threshold = 0.005;
        int ppm = 500;

        updateDatabaseWithMergedSpectra(databaseFileName, ppm);
    }

    /**
     * Generates a combined spectrum for each peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  databaseFileName database to pull data from
     * @param  ppm parts per million value
     */
    private static void updateDatabaseWithMergedSpectra(String databaseFileName, int ppm) throws IOException {

        int numOfPeptideIds = 16474;

        spectrumDbUtil copyDB = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");
        spectrumDbUtil specDB = new spectrumDbUtil("testLibDuplicateSpectra.db");
        specDB.connect();
        copyDB.connect();

        copyDB.createNewMergedSpectraTable();

        for(int peptideID = 0; peptideID < numOfPeptideIds; peptideID++){

            fileManager fm = new fileManager();
            spectrumMerger sm =  new spectrumMerger(peptideID);

            PeakList[] individualPeptidePeakLists = specDB.getSpectrumPeakIntensityByProteinID2(peptideID);

            sm.addPeaksToSpecArray(individualPeptidePeakLists);
            sm.mergePeaksByThreshold(ppm);

            // Add the spec array to the database
            copyDB.addMergedPeptideRow(peptideID, sm.getSpecArr(), sm.getNumberOfNonZeroElements());

            //System.out.println("Done! Peptide ID: " + peptideID);
        }
        specDB.disconnect();
        copyDB.disconnect();



    }

    /**
     * Generates a combined spectrum for each peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  databaseFileName database to pull data from
     * @param  ppm parts per million value
     */
    public void updateDatabaseWithMergedSpectraMultiThread(String databaseFileName, int ppm, int beginPeptide, int endPeptide) throws IOException {

        System.out.println("begin: " + beginPeptide);
        for(int peptideID = beginPeptide; peptideID < endPeptide; peptideID++){

            fileManager fm = new fileManager();
            spectrumMerger sm =  new spectrumMerger(peptideID);

            PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

            sm.addPeaksToSpecArray(individualPeptidePeakLists);
            sm.mergePeaksByThreshold(ppm);

            // Add the spec array to the database
            if(copyDBMT.isConnected()) {
                copyDBMT.addMergedPeptideRow(peptideID, sm.getSpecArr(), sm.getNumberOfNonZeroElements());
            }else{
                copyDBMT.connect();
                copyDBMT.addMergedPeptideRow(peptideID, sm.getSpecArr(), sm.getNumberOfNonZeroElements());
            }

            System.out.println("Done! Peptide ID: " + peptideID);
        }

    }

    /**
     * Generates a combined spectrum for each peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  databaseFileName database to pull data from
     * @param  ppm parts per million value
     */
    private static void generateAllPeptideSpectrumMerged(String databaseFileName, int ppm) throws IOException {

        int numOfPeptideIds = 10;

        String dirName = "outputCSV/MergedSpectrum";

        for(int peptideID = 0; peptideID < numOfPeptideIds; peptideID++){

            fileManager fm = new fileManager();
            spectrumMerger sm =  new spectrumMerger(peptideID);

            PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

            sm.addPeaksToSpecArray(individualPeptidePeakLists);
            sm.mergePeaksByThreshold(ppm);
            fm.writeSpectrumArrayToCSV(sm.getSpecArr(), peptideID, dirName);

            System.out.println("Done");
        }

    }

    /**
     * Generates a combined spectrum for each peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  databaseFileName database to pull data from     */
    private static void generateAllPeptideSpectrumCombined(String databaseFileName) throws IOException {

        int numOfPeptideIds = 10;

        //DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        //Date date = new Date();

        String dirName = "outputCSV/CombinedSpectrum";

        for(int peptideID = 0; peptideID < numOfPeptideIds; peptideID++){

            fileManager fm = new fileManager();
            spectrumMerger sm =  new spectrumMerger(peptideID);

            PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

            sm.addPeaksToSpecArray(individualPeptidePeakLists);
            fm.writeSpectrumArrayToCSV(sm.getSpecArr(), peptideID, dirName);

            System.out.println("Done");
        }

    }

    /**
     * Generates a combined spectrum for the given peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  peptideID the ID of the peptide to generate a combined spectrum
     * @param  databaseFileName database to pull data from     */
    public static void generateSinglePeptideSpectrumMerged(int peptideID, String databaseFileName, int ppm) throws IOException {

        fileManager fm = new fileManager();
        spectrumMerger sm =  new spectrumMerger(peptideID);

        //DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        //Date date = new Date();
        String dirName = "outputCSV";

        PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

        sm.addPeaksToSpecArray(individualPeptidePeakLists);
        sm.mergePeaksByThreshold(ppm);

        fm.writeSpectrumArrayToCSV(sm.getSpecArr(), peptideID, dirName);

        System.out.println("Done");

    }

    /**
     * Generates a combined spectrum for the given peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  peptideID the ID of the peptide to generate a combined spectrum
     * @param  databaseFileName database to pull data from     */
    public static void generateSinglePeptideSpectrumCombined(int peptideID, String databaseFileName) throws IOException {
        fileManager fm = new fileManager();
        spectrumMerger sm =  new spectrumMerger(peptideID);
        String dirName = "outputCSV";

        PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

        sm.addPeaksToSpecArray(individualPeptidePeakLists);

        fm.writeSpectrumArrayToCSV(sm.getSpecArr(), peptideID, dirName);
        System.out.println("Done");
    }

    /**
     * Generates the individual spectrum's for the given peptide in the database and outputs the spectrum as a csv
     * file labeled: individualSpectrum[PeptideID].txt
     *
     * @param  peptideID the ID of the peptide to generate a combined spectrum
     * @param  databaseFileName database to pull data from
     */
    public static void generateSingleIndividualCSV(int peptideID, String databaseFileName) throws IOException {
        fileManager fm = new fileManager();
        spectrumMerger sm =  new spectrumMerger(peptideID);
        String dirName = "outputCSV/IndividualSpectrum";

        PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

        //Generate sample individual csv's for each uncombined peak in Peaklist
        fm.writePeakListToCSV(individualPeptidePeakLists, individualPeptidePeakLists.length, dirName);
        System.out.println("Done");
    }


    private static PeakList[] getPeakListFromDatabase(int peptideID, String databaseFileName){
        spectrumDbUtil specDB = new spectrumDbUtil(databaseFileName);

        PeakList[] individualPeptidePeakLists;
        if(specDB.connect()) {
            individualPeptidePeakLists = specDB.getSpectrumPeakIntensityByProteinID2(peptideID);
            specDB.disconnect();
        }else{
            System.out.println("Error Opening connection to Database");
            return null;
        }

        if(individualPeptidePeakLists == null){
            System.out.println("individualPeptidePeakLists is null: exiting");
            System.exit(1);
        }

        return individualPeptidePeakLists;
    }

}
