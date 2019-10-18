package SpectrumMerger;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import javafx.util.Pair;

public class spectrumDbUtil {

    String databaseFileName = "testLibDuplicateSpectra.db";
    Connection specDbConnection;

    public spectrumDbUtil(String databaseFileName){
        this.databaseFileName = databaseFileName;
    }

    public spectrumDbUtil(){}


    public PeakList getSpectrumPeakIntensityByProteinID(int proteinID){

        FloatBuffer floatBuf;
        float[] peakMzArray;
        float[] peakIntensityArray;

        PeakList list = new PeakList();

        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT * FROM SpectraTable WHERE peptideID = " + proteinID)){

            while (rs.next()) {

                byte[] peakMzByteArr = rs.getBytes("peakMZ");
                byte[] intensityByteArr = rs.getBytes("peakIntensity");

                //convert PeakMz byteArray to float array
                floatBuf = ByteBuffer.wrap(peakMzByteArr).order(ByteOrder.BIG_ENDIAN).asFloatBuffer();
                peakMzArray = new float[floatBuf.remaining()];
                floatBuf.get(peakMzArray);

                //convert PeakIntensity byteArray to float array
                floatBuf = ByteBuffer.wrap(intensityByteArr).order(ByteOrder.BIG_ENDIAN).asFloatBuffer();
                peakIntensityArray = new float[floatBuf.remaining()];
                floatBuf.get(peakIntensityArray);

                //combine 2 arrays into peaklist array
                for(int i=0;i<peakMzArray.length;i++){
                    Peak p = new Peak((double)peakMzArray[i], (double)peakIntensityArray[i]);
                    list.addPeak(p);
                }
            }
            return list;
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
        return null;
    }
    public PeakList[] getSpectrumPeakIntensityByProteinID2(int proteinID){

        FloatBuffer floatBuf;
        float[] peakMzArray;
        float[] peakIntensityArray;

        PeakList[] listOfPeaks;
        int numberOfDBHits = 0;

        //get size of database query
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT COUNT(*) FROM SpectraTable WHERE peptideID = " + proteinID);
        ){
            numberOfDBHits = rs.getInt(1);
        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return null;
        }

        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT * FROM SpectraTable WHERE peptideID = " + proteinID);
        ){

            listOfPeaks = new PeakList[numberOfDBHits];

            int counter = 0;
            while (rs.next()) {

                byte[] peakMzByteArr = rs.getBytes("peakMZ");
                byte[] intensityByteArr = rs.getBytes("peakIntensity");

                //convert PeakMz byteArray to float array
                floatBuf = ByteBuffer.wrap(peakMzByteArr).order(ByteOrder.BIG_ENDIAN).asFloatBuffer();
                peakMzArray = new float[floatBuf.remaining()];
                floatBuf.get(peakMzArray);

                //convert PeakIntensity byteArray to float array
                floatBuf = ByteBuffer.wrap(intensityByteArr).order(ByteOrder.BIG_ENDIAN).asFloatBuffer();
                peakIntensityArray = new float[floatBuf.remaining()];
                floatBuf.get(peakIntensityArray);

                //combine 2 arrays into array of pairs
                listOfPeaks[counter] = new PeakList();
                for(int i=0;i<peakMzArray.length;i++){
                    Peak p = new Peak((double)peakMzArray[i], (double)peakIntensityArray[i]);

                    listOfPeaks[counter].addPeak(p);

                }
                counter++;
            }
            return listOfPeaks;
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
        return null;
    }

    public int[] getAllPossibleMasses(int proteinID) {
        //get size of database query
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT COUNT(*) FROM Peptidetable WHERE peptideID = " + proteinID);
        ){

        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return null;
        }

        return null;
    }


    public ResultSet execute(String query){
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery(query)){

            while (rs.next()) {
                //System.out.println(rs.getInt(1));
            }

            return rs;

        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
        return null;
    }


    public void getTableNames() throws SQLException {
        ResultSet rs = specDbConnection.getMetaData().getTables(null, null, null, null);
        while (rs.next()) {
            System.out.println(rs.getString("TABLE_NAME"));
        }
        System.out.println(specDbConnection.getSchema());

    }

    public boolean connect() {
        try {
            // db parameters
            String url = "jdbc:sqlite:" + databaseFileName;
            // create a connection to the database
            specDbConnection = DriverManager.getConnection(url);

            System.out.println("Connection to SQLite has been established.");

        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return false;
        }
        return true;
    }

    public boolean disconnect() {
        try {
            if (specDbConnection != null) {
                specDbConnection.close();
                System.out.println("Connection to SQLite has been closed.");
            }
        } catch (SQLException ex) {
            System.out.println(ex.getMessage());
            return false;
        }

        return true;
    }

    public double[] getPeptideSpectraArray(){
        if(specDbConnection == null) return null;

        // Query db for spectraTable
        // connvert
        return null;
    }

    public static void main(String[] args) throws IOException, SQLException {
        spectrumDbUtil specDB = new spectrumDbUtil();
        if(specDB.connect() == true) {
            PeakList[] list = specDB.getSpectrumPeakIntensityByProteinID2(1);
            System.out.println(list.length);
            specDB.disconnect();
        }
    }
}
