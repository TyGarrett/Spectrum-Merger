package SpectrumMerger;

import javafx.util.Pair;
import org.sqlite.SQLiteConfig;

import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.sql.*;
import java.util.LinkedList;

/**
 * Useful functions to handle database for merging spectra
 */
public class SpectrumMergerDataBaseUtil {

    private String databaseFileName;
    private Connection specDbConnection;

    private String MergedSpectrumTable = "MergedSpectraTable";
    private String SpectraTable = "SpectraTable";

    public SpectrumMergerDataBaseUtil(String dbFileName){
        this.databaseFileName = dbFileName;
    }

    /**
     * Updates MergedSpectraTable to SpectraTable
     */
    void changeMergedTableNameToSpectrumTable(){
        try (Statement stmt  = specDbConnection.createStatement()
        ){
            int deletedRows = stmt.executeUpdate("ALTER TABLE "+MergedSpectrumTable+" RENAME TO "+SpectraTable);
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Deletes SpectraTable
     */
    void deleteSpectraTable(){
        try (Statement stmt  = specDbConnection.createStatement()
        ){
            int deletedRows = stmt.executeUpdate("DROP TABLE IF EXISTS "+SpectraTable);

        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Creates new MergedSpectraTable and deletes existing one
     */
    void createNewMergedSpectraTable(){

        // SQL statement for creating a new table
        String sql = "CREATE TABLE IF NOT EXISTS "+MergedSpectrumTable+" (\n"
                + "    peptideID integer PRIMARY KEY,\n"
                + "    peakMZ Binary,\n"
                + "    peakIntensity Binary,\n"
                + "    massKey integer\n"
                + ");";

        String sqlDrop = "DROP TABLE IF EXISTS "+MergedSpectrumTable;

        try (Statement stmt  = specDbConnection.createStatement()
        ){
            stmt.execute(sqlDrop);
            stmt.execute(sql);
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Insert a a batch of merged peptide spectra into the database
     */
    public synchronized void insertBatch(int startPeptideID, LinkedList<float[]> listOfMzArrays,
                                         LinkedList<float[]>listOfIntensityArrays,  LinkedList<Integer> massKeys)
            throws SQLException {

        String SQL = "INSERT INTO "+MergedSpectrumTable+"(peptideID,peakMz,peakIntensity,massKey) VALUES(?,?,?,?)";
        PreparedStatement pstmt = specDbConnection.prepareStatement(SQL);
        specDbConnection.setAutoCommit(false);

        for (int i = 0; i < listOfMzArrays.size(); i++){
            byte[] peakMzByteArr;
            byte[] peakIntensityByteArr;

            peakMzByteArr = FloatArray2ByteArray(listOfMzArrays.get(i));
            peakIntensityByteArr = FloatArray2ByteArray(listOfIntensityArrays.get(i));

            int peptideID = startPeptideID + i;

            pstmt.setInt(1, peptideID);
            pstmt.setBytes(2, peakMzByteArr);
            pstmt.setBytes(3, peakIntensityByteArr);
            pstmt.setInt(4, massKeys.get(i));

            // Add above SQL statement in the batch.
            pstmt.addBatch();
        }
        //Create an int[] to hold returned values
        int[] count = pstmt.executeBatch();

        //Explicitly commit statements to apply changes
        specDbConnection.commit();

    }

    /**
     * returns the total number of unique peptides in the SpectraTable
     */
    int getNumberOfUniquePeptideIDsFromSpectraTable(){
        int numPeptides = 0;
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT * FROM SpectraTable ORDER BY peptideID DESC LIMIT 0, 1")
        ){
            numPeptides = rs.getInt(1);
        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return 0;
        }
        return numPeptides;
    }

    /**
     * @param proteinID the unique id of the protein spectra being merged
     * @return total number of unique peptides in the SpectraTable
     */
    int getNumberOfSpectraWithProteinID(int proteinID){
        int numberOfDBHits;

        //get size of database query
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT COUNT(*) FROM SpectraTable WHERE peptideID = " + proteinID)
        ){
            numberOfDBHits = rs.getInt(1);
        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return -1;
        }
        return numberOfDBHits;
    }

    /**
     * generates the spectrumMerger object with all spectra with peptideID
     *
     * @param proteinID the unique id of the protein spectra being merged
     * @return spectrumMerger object with all spectra with proteinID included
     */
    spectrumMerger generateSpectrumMerger(int proteinID){

        FloatBuffer floatBuf;
        float[] peakMzArray;
        float[] peakIntensityArray;

        spectrumMerger sm = new spectrumMerger(proteinID);
        int numberOfDBHits = getNumberOfSpectraWithProteinID(proteinID);

        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT * FROM SpectraTable WHERE peptideID = " + proteinID)
        ){

            while (rs.next()) {

                LinkedList<Pair<Float, Float>> pairList =  new LinkedList<Pair<Float, Float>>();
                float maxIntensity = 0;

                byte[] peakMzByteArr = rs.getBytes("peakMZ");
                byte[] intensityByteArr = rs.getBytes("peakIntensity");

                sm.setMassKey(rs.getInt("massKey"));

                //convert PeakMz byteArray to float array
                floatBuf = ByteBuffer.wrap(peakMzByteArr).order(ByteOrder.BIG_ENDIAN).asFloatBuffer();
                peakMzArray = new float[floatBuf.remaining()];
                floatBuf.get(peakMzArray);

                //convert PeakIntensity byteArray to float array
                floatBuf = ByteBuffer.wrap(intensityByteArr).order(ByteOrder.BIG_ENDIAN).asFloatBuffer();
                peakIntensityArray = new float[floatBuf.remaining()];
                floatBuf.get(peakIntensityArray);

                for(int i=0;i<peakIntensityArray.length;i++){
                    if(peakIntensityArray[i] > maxIntensity){
                        maxIntensity = peakIntensityArray[i];
                    }
                }
                for(int i=0;i<peakIntensityArray.length;i++){
                    pairList.add(new Pair(peakMzArray[i] ,peakIntensityArray[i]/maxIntensity));
                }

                sm.addPairToSpecArray(pairList);

            }
            return sm;
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
        return null;
    }

    /**
     * connects to the database with properties to boost insert performance
     *
     * @return true if database is connected, false otherwise
     */
    boolean connect() {
        try {
            String url = "jdbc:sqlite:" + databaseFileName;
            final int cacheSize = 100000 / 6;
            final int pageSize = 4096;
            SQLiteConfig config = new SQLiteConfig();
            //optimize for multiple connections that can share data structures
            config.setSharedCache(true);
            config.setCacheSize(cacheSize);
            config.setPageSize(pageSize);
            config.setJournalMode(SQLiteConfig.JournalMode.TRUNCATE);
            config.enableFullSync(false);
            config.enableRecursiveTriggers(false);
            config.setLockingMode(SQLiteConfig.LockingMode.NORMAL);
            config.setSynchronous(SQLiteConfig.SynchronousMode.OFF); //TODO may be dangerous on some systems to have off
            specDbConnection = DriverManager.getConnection(url,config.toProperties());
        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return false;
        }
        return true;
    }

    /**
     * disconnects from the database
     */
    void disconnect() {
        try {
            if (specDbConnection != null) {
                specDbConnection.close();
            }
        } catch (SQLException ex) {
            System.out.println(ex.getMessage());
        }
    }

    /**
     * converts a float array into a binary array for the database
     *
     * @param floatArray the array to be converted
     * @return the converted byte array
     */
    public static byte[] FloatArray2ByteArray(float[] floatArray) {

        ByteArrayOutputStream bous = new ByteArrayOutputStream();
        DataOutputStream dout = new DataOutputStream(bous);
        try {
            for (float d : floatArray) {
                dout.writeFloat(d);
                //bous.write(Float.floatToIntBits(d));
            }
            dout.close();
        } catch(IOException e){
            System.out.println(e);
        }
        return bous.toByteArray();
    }
}