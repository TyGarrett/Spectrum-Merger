package SpectrumMerger;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import javafx.util.Pair;
import org.apache.commons.math3.complex.Complex;
import org.jboss.util.Null;

import java.io.FileWriter;
import java.io.IOException;

public class fileManager {

    public fileManager(){

    }
    public void writeFFTtoCSV(Complex[] comp) throws IOException {
        FileWriter csvWriter = new FileWriter("/home/ty/jupyterNotebooks/FFTSpectrum.csv");
        csvWriter.append("real");
        csvWriter.append(",");
        csvWriter.append("Imaginary");
        csvWriter.append(",");
        csvWriter.append("combined");
        csvWriter.append("\n");

        int counter = 0;
        for (Complex c : comp) {
            double dist = Math.pow(c.getReal(),2) + Math.pow(c.getImaginary(), 2);
            dist = Math.sqrt(dist);

            csvWriter.append(round(c.getReal(),5)+","+round(c.getImaginary(),5)+","+round(dist,5));
            csvWriter.append("\n");

            counter++;
        }

        csvWriter.flush();
        csvWriter.close();
    }

    public void writeCombinedSpectrumtoCSV(int[] specArr, int peptideID) throws IOException {
        FileWriter csvWriter = new FileWriter("/home/ty/jupyterNotebooks/CombinedSpectrum"+peptideID+".csv");
        int counter = 0;
        for (int intensity : specArr) {
            csvWriter.append(intensity+",");
            csvWriter.append("\n");

            counter++;
        }

        csvWriter.flush();
        csvWriter.close();

        System.out.println("Wrote file{}".format(String.valueOf(peptideID)));
    }

    public void writeIndividualSpectrumtoCSV(PeakList[] list, int numFilesToCreate) throws IOException {

        int counter = 0;
        for (PeakList pl : list) {

            if(counter == numFilesToCreate) return;

            FileWriter csvWriter = new FileWriter("/home/ty/jupyterNotebooks/Spec/IndividualSpectrum" + counter + ".csv");
            csvWriter.append("Mass");
            csvWriter.append(",");
            csvWriter.append("Intensity");
            csvWriter.append("\n");

            for (Peak p : pl.getSortedPeaks(false)) {

                double intensity = p.getIntensity();
                double M2z = p.getM2z();

                csvWriter.append(M2z+","+intensity);
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
