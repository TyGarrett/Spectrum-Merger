package SpectrumMerger;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;

import java.io.FileWriter;
import java.io.IOException;

public class fileManager {

    fileManager(){}

    void writeSpectrumArrayToCSV(float[] specArr, int peptideID, String dirPath) throws IOException {

        FileWriter csvWriter = new FileWriter(dirPath+"/Spectrum["+peptideID+"].csv");

        for (float intensity : specArr) {
            csvWriter.append(String.valueOf(intensity)).append(",");
            csvWriter.append("\n");
        }

        csvWriter.flush();
        csvWriter.close();

        System.out.println(peptideID);
    }

    void writePeakListToCSV(PeakList[] list, int numFilesToCreate, String dirPath) throws IOException {

        int counter = 0;
        for (PeakList pl : list) {

            if(counter == numFilesToCreate) return;

            FileWriter csvWriter = new FileWriter(dirPath+"/IndividualSpectrum["+counter+"].csv");

            csvWriter.append("Mass");
            csvWriter.append(",");
            csvWriter.append("Intensity");
            csvWriter.append("\n");

            for (Peak p : pl.getSortedPeaks(false)) {

                double intensity = p.getIntensity();
                double M2z = p.getM2z();

                csvWriter.append(String.valueOf(M2z)).append(",").append(String.valueOf(intensity));
                csvWriter.append("\n");
            }

            csvWriter.flush();
            csvWriter.close();

            counter++;
        }
    }

    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        long factor = (long) Math.pow(10, places);
        value = value * factor;
        long tmp = Math.round(value);
        return (double) tmp / factor;
    }
}
