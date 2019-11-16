package SpectrumMerger;

import javafx.util.Pair;
import java.sql.SQLException;
import java.util.LinkedList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Independent threads to handle merging spectra and adding to database
 */
class SpectrumMergeRunnable implements Runnable {

    private static final double SECOND_CONVERSION = 1000000000.0;

    int threadCount;
    SpectrumMergerDataBaseUtil multiDB;
    int numberOfPeptidesPerThread;
    int totalProteins;
    int peptideIdStartIndex;
    int peptideIdEndIndex;

    int ppm = 500;
    int batchSize = 500;

    /**
     * Generates the individual spectrum's for the given peptide in the database and outputs the spectrum as a csv
     * file labeled: individualSpectrum[PeptideID].txt
     *
     * @param  threadCount 0 - number of threads
     * @param  multiDB open database file
     * @param  numberOfPeptidesPerThread number of peptides to merge
     * @param  totalProteins total number of peptides in database SpectraTable
     */
    SpectrumMergeRunnable(int threadCount, SpectrumMergerDataBaseUtil multiDB, int numberOfPeptidesPerThread,
                          int totalProteins){

        this.threadCount = threadCount;
        this.multiDB = multiDB;
        this.numberOfPeptidesPerThread = numberOfPeptidesPerThread;
        this.totalProteins = totalProteins;

        //calculates the start and end Peptide ID. start: inclusive, end exclusive
        this.peptideIdStartIndex = threadCount*numberOfPeptidesPerThread;
        this.peptideIdEndIndex = ++threadCount*numberOfPeptidesPerThread;
    }

    public void run() {
        try {
            MergeSpectra();
        } catch (Exception e) {
            System.out.println("Exception is caught");
            System.out.println(e);
        }
    }

    /**
     * Generates the individual spectrum's for the given peptide in the database and outputs the spectrum as a csv
     * file labeled: individualSpectrum[PeptideID].txt
     */
    private void MergeSpectra() throws SQLException {
        // Initialize the empty list for each batch of database entries
        LinkedList<float[]> listOfMzArrays = new LinkedList<>();
        LinkedList<float[]> listOfIntensityArrays = new LinkedList<>();
        LinkedList<Integer> listOfMassKeys = new LinkedList<>();

        long startTime = System.nanoTime();
        int peptideIDAtStartOfBatch = peptideIdStartIndex;
        for(int peptideID = peptideIdStartIndex; peptideID < peptideIdEndIndex; peptideID++){

            spectrumMerger sm = multiDB.generateSpectrumMerger(peptideID); // add peaks together into one spectra
            sm.mergePeaksByThreshold(ppm); //merge peaks together according to ppm

            Pair<float[], float[]> massIntensityPair = sm.getMassIntensityArrays();

            // populate the lists with appropriate database values
            listOfMzArrays.add(massIntensityPair.getKey());
            listOfIntensityArrays.add(massIntensityPair.getValue());
            listOfMassKeys.add(sm.getMassKey());

            // once list size = batch size then write to database
            if(listOfMzArrays.size() >= batchSize || peptideID >= totalProteins || peptideID == peptideIdEndIndex - 1){

                multiDB.insertBatch(peptideIDAtStartOfBatch, listOfMzArrays, listOfIntensityArrays, listOfMassKeys);

                // reinitialize the empty lists for each batch of database entries
                listOfMzArrays = new LinkedList<>();
                listOfIntensityArrays = new LinkedList<>();
                listOfMassKeys = new LinkedList<>();

                peptideIDAtStartOfBatch = peptideID + 1;

                // break if the thread has reached end if
                if(peptideID >= totalProteins){
                    break;
                }
            }
        }
        int threadDuration = (int)((System.nanoTime() - startTime)/SECOND_CONVERSION);
        System.out.println("thread " +threadCount+" total time (s): " + threadDuration);
        System.out.println("    ProteinIDs: " + peptideIdStartIndex + " - " + peptideIdEndIndex);
    }
}

/**
 * Creates database and handles threads
 */
class SpectrumMergeThreadHandler {
    public static void main(String[] args) throws InterruptedException {

        // create database connection to read and write
        SpectrumMergerDataBaseUtil dataBase = new SpectrumMergerDataBaseUtil("lib_2908redundant.db");
        dataBase.connect();

        // add merged table to spectra table, remove existing mergedTable if exists
        dataBase.createNewMergedSpectraTable();

        ExecutorService es = Executors.newCachedThreadPool();
        int numberOfThreads = 20;
        int greatestPeptideID = dataBase.getNumberOfUniquePeptideIDsFromSpectraTable();
        System.out.println("Total number of peptides: " + (greatestPeptideID + 1));

        // handles decimal point (upper bound) from dividing total peptides by thread count
        int numberOfPeptidesPerThread = greatestPeptideID/numberOfThreads + 1;

        // create threads to handle merging of spectra
        for(int i=0;i<numberOfThreads;i++) {
            Thread object = new Thread(new SpectrumMergeRunnable(i, dataBase, numberOfPeptidesPerThread, greatestPeptideID));
            es.execute(object);
        }

        es.shutdown();
        boolean finished = es.awaitTermination(10, TimeUnit.MINUTES);

        // when all threads are finished remove original SpectraTable and replace with Mere
        if(finished){
            //TODO: DOESN'T WORK... IDK WHY????
            System.out.println("table alterations");
            dataBase.deleteSpectraTable();
            dataBase.changeMergedTableNameToSpectrumTable();
        }
        dataBase.disconnect();
    }
}



