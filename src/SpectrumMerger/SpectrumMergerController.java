package SpectrumMerger;

import edu.scripps.pms.util.spectrum.PeakList;

import java.io.IOException;
import java.sql.SQLException;
import java.util.LinkedList;

class SpectrumMergerHandler {
    spectrumDbUtil copyDBMT = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");

    public SpectrumMergerHandler(){

    }

    public static void main(String[] args) throws IOException, SQLException {
        String databaseFileName = "testLibDuplicateSpectra.db";
        double threshold = 0.005;
        int ppm = 500;

        timeComplexityOfBatch(databaseFileName, ppm);
    }

    /**
     * Generates a combined spectrum for each peptide in the database and outputs the value as a csv file labeled:
     * combinedSpectrum[PeptideID].txt
     *
     * @param  databaseFileName database to pull data from
     * @param  ppm parts per million value
     */
    private static void updateDatabaseWithMergedSpectraBatch(String databaseFileName, int ppm) throws SQLException {
        int numOfPeptideIds = 16474;
        int batchSize = 30;

        LinkedList<float[]> listOfSpecArrays = new LinkedList<float[]>();
        LinkedList<Integer> numElems = new LinkedList<Integer>();

        spectrumDbUtil copyDB = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");
        copyDB.connect();
        copyDB.createNewMergedSpectraTable();

        long startTime = System.nanoTime();

        int startPeptideID = 0;
        for(int peptideID = 0; peptideID < numOfPeptideIds; peptideID++){

            spectrumMerger sm =  new spectrumMerger(peptideID);

            PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);

            sm.addPeaksToSpecArray(individualPeptidePeakLists);
            sm.mergePeaksByThreshold(ppm);/*
 */
            listOfSpecArrays.add(sm.getSpecArr());
            numElems.add(sm.getNumberOfNonZeroElements());

            if(listOfSpecArrays.size() >= batchSize || peptideID >= numOfPeptideIds ){
                //copyDB.insertBatch(startPeptideID, listOfSpecArrays, numElems);

                listOfSpecArrays = new LinkedList<float[]>();
                numElems = new LinkedList<Integer>();
                startPeptideID = peptideID + 1;

                long endTime = System.nanoTime();
                long duration = (endTime - startTime);
                System.out.println("Sec's left: " + (int)((duration/1000000000.0)/(((double)peptideID)/numOfPeptideIds) - (duration/1000000000.0)));
                System.out.println("percent: " + ((double)peptideID)/numOfPeptideIds*100);
            }
        }
        copyDB.disconnect();
    }

    private static void timeComplexityOfBatch(String databaseFileName, int ppm) throws SQLException {
        long durationSM = 0;
        long durationGP = 0;
        long durationAP = 0;
        long durationMP = 0;
        long durationDB = 0;

        //int numOfPeptideIds = 16474;
        int numOfPeptideIds = 100;
        int batchSize = 30;

        LinkedList<float[]> listOfSpecArrays = new LinkedList<float[]>();
        LinkedList<Integer> numElems = new LinkedList<Integer>();

        spectrumDbUtil copyDB = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");
        copyDB.connect();
        copyDB.createNewMergedSpectraTable();

        long startTime = System.nanoTime();

        int startPeptideID = 0;
        for(int peptideID = 0; peptideID < numOfPeptideIds; peptideID++){

            long startSM = System.nanoTime();
            //spectrumMerger sm =  new spectrumMerger(peptideID);
            spectrumMerger sm = getSpectrumMergerForPeptide(peptideID, databaseFileName);
            durationSM += (System.nanoTime() - startSM);

            long start2 = System.nanoTime();
            //PeakList[] individualPeptidePeakLists = getPeakListFromDatabase(peptideID, databaseFileName);
            durationGP += (System.nanoTime() - start2);

            long start3 = System.nanoTime();
            //sm.addPeaksToSpecArray(individualPeptidePeakLists);
            durationAP += (System.nanoTime() - start3);

            long start4 = System.nanoTime();
            sm.mergePeaksByThreshold(ppm);
            durationMP += (System.nanoTime() - start4);

            listOfSpecArrays.add(sm.getSpecArr());
            numElems.add(sm.getNumberOfNonZeroElements());

            if(listOfSpecArrays.size() >= batchSize || peptideID >= numOfPeptideIds ){
                long start5 = System.nanoTime();
                copyDB.insertBatch(startPeptideID, listOfSpecArrays, numElems);
                durationDB += (System.nanoTime() - start5);

                listOfSpecArrays = new LinkedList<float[]>();
                numElems = new LinkedList<Integer>();
                startPeptideID = peptideID + 1;

                long endTime = System.nanoTime();
                long duration = (endTime - startTime);
                System.out.println("Sec's left: " + (int)((duration/1000000000.0)/(((double)peptideID)/numOfPeptideIds) - (duration/1000000000.0)));
                System.out.println("percent: " + ((double)peptideID)/numOfPeptideIds*100);
            }
        }
        System.out.println("duration SpecMerge Construct: " + (int)((durationSM/1000000.0)));
        System.out.println("duration peaklist from DB   : " + (int)((durationGP/1000000.0)));
        System.out.println("duration add peaks to SM    : " + (int)((durationAP/1000000.0)));
        System.out.println("duration merge peaks in SM  : " + (int)((durationMP/1000000.0)));
        System.out.println("write peaks batch to DB     : " + (int)((durationDB/1000000.0)));
        System.out.println("");

        System.out.println("total time : " + (int)((System.nanoTime() - startTime)/1000000.0));

        copyDB.disconnect();
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
    private static spectrumMerger getSpectrumMergerForPeptide(int peptideID, String databaseFileName){
        spectrumDbUtil specDB = new spectrumDbUtil(databaseFileName);

        spectrumMerger sm;

        if(specDB.connect()) {
            sm = specDB.generateSpectrumMerger(peptideID);
            specDB.disconnect();
        }else{
            System.out.println("Error Opening connection to Database");
            return null;
        }

        if(sm == null){
            System.out.println("individualPeptidePeakLists is null: exiting");
            System.exit(1);
        }

        return sm;
    }


    public void updateDatabaseWithMergedSpectraMultiThread(String databaseFileName, int ppm, int start, int end, SpectrumMergerDataBaseUtil multiDB, int batchSize, int numOfPeptideIds) throws SQLException {

        LinkedList<float[]> listOfSpecArrays = new LinkedList<float[]>();
        LinkedList<Integer> numElems = new LinkedList<Integer>();

        long startTime = System.nanoTime();
        int startPeptideID = start;
        for(int peptideID = start; peptideID < end; peptideID++){
            spectrumMerger sm = getSpectrumMergerForPeptide(peptideID, databaseFileName);

            sm.mergePeaksByThreshold(ppm);

            listOfSpecArrays.add(sm.getSpecArr());
            numElems.add(sm.getNumberOfNonZeroElements());

            if(listOfSpecArrays.size() >= batchSize || peptideID >= numOfPeptideIds ){
                //multiDB.insertBatch(peptideID, listOfSpecArrays, numElems);

                listOfSpecArrays = new LinkedList<float[]>();
                numElems = new LinkedList<Integer>();
                startPeptideID = peptideID + 1;

            }
        }
        System.out.println("total time (" + start+") : " + (int)((System.nanoTime() - startTime)/1000000.0));
    }
}
