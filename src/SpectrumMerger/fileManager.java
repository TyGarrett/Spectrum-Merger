package SpectrumMerger;

import org.apache.commons.math3.complex.Complex;

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

    public void writeCombinedSpectrumtoCSV(double[] comp) throws IOException {
        FileWriter csvWriter = new FileWriter("/home/ty/jupyterNotebooks/CombinedSpectrum.csv");

        int counter = 0;
        for (double c : comp) {
            csvWriter.append(round(c,5)+",");
            csvWriter.append("\n");

            counter++;
        }

        csvWriter.flush();
        csvWriter.close();
    }

    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        long factor = (long) Math.pow(10, places);
        value = value * factor;
        long tmp = Math.round(value);
        return (double) tmp / factor;
    }
}
