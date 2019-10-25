package SpectrumMerger;


// Java code for thread creation by implementing
// the Runnable Interface
class MultiThreadSpectrumMerger implements Runnable {

    int count = 0;
    MultiThreadSpectrumMerger(int c){
        count = c;
    }

    public void run() {
        try {
            // Displaying the thread that is running
            System.out.println("Thread " +
                    Thread.currentThread().getId() +
                    " is running");

            String databaseFileName = "testLibDuplicateSpectra.db";
            double threshold = 0.005;
            int ppm = 500;

            SpectrumMergerHandler smh = new SpectrumMergerHandler();
            System.out.println(count);
            smh.updateDatabaseWithMergedSpectraMultiThread(databaseFileName, ppm, count*1000, ++count*1000);

        } catch (Exception e) {
            // Throwing an exception
            System.out.println("Exception is caught");
        }
    }
}

// Main Class
class Multithread {
    public static void main(String[] args) {
        spectrumDbUtil copyDB = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");
        copyDB.connect();
        copyDB.createNewMergedSpectraTable();
        copyDB.disconnect();
        System.out.println("DROPPING!!!!!!!!!!!!!!");

        int n = 8; // Number of threads
        for (int i = 0; i < 8; i++) {
            Thread object = new Thread(new MultiThreadSpectrumMerger(i));
            object.start();
        }
    }
}


