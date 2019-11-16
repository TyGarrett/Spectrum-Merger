package SpectrumMerger;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.sql.*;
import java.util.LinkedList;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import javafx.util.Pair;
import org.sqlite.SQLiteConfig;

public class spectrumDbUtil {

    private String databaseFileName;
    private Connection specDbConnection;
    private static boolean isConnected = false;

    public spectrumDbUtil(String dbFileName){
        this.databaseFileName = dbFileName;
    }

    public void changeTableName(String tableToChange, String name){
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs   = stmt.executeQuery("ALTER TABLE "+tableToChange+" RENAME TO "+name)
        ){

        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }

    public void deleteTable(String tableToDelete){
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs   = stmt.executeQuery("DROP TABLE"+tableToDelete);
        ){

        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }

    PeakList[] getSpectrumPeakIntensityByProteinID2(int proteinID){

        FloatBuffer floatBuf;
        float[] peakMzArray;
        float[] peakIntensityArray;

        PeakList[] listOfPeaks;
        int numberOfDBHits;

        //get size of database query
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT COUNT(*) FROM SpectraTable WHERE peptideID = " + proteinID)
        ){
            numberOfDBHits = rs.getInt(1);
        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return null;
        }

        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT * FROM SpectraTable WHERE peptideID = " + proteinID)
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
                    Peak p = new Peak(peakMzArray[i], peakIntensityArray[i]);

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



    synchronized spectrumMerger generateSpectrumMerger(int proteinID){

        FloatBuffer floatBuf;
        float[] peakMzArray;
        float[] peakIntensityArray;

        spectrumMerger sm = new spectrumMerger(proteinID);

        int numberOfDBHits;

        //get size of database query
        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT COUNT(*) FROM SpectraTable WHERE peptideID = " + proteinID)
        ){
            numberOfDBHits = rs.getInt(1);
        } catch (SQLException e) {
            System.out.println(e.getMessage());
            return null;
        }

        try (Statement stmt  = specDbConnection.createStatement();
             ResultSet rs    = stmt.executeQuery("SELECT * FROM SpectraTable WHERE peptideID = " + proteinID)
        ){

            boolean setMassKey = false;
            while (rs.next()) {

                LinkedList<Pair<Float, Float>> pairList =  new LinkedList<Pair<Float, Float>>();
                float maxIntensity = 0;

                byte[] peakMzByteArr = rs.getBytes("peakMZ");
                byte[] intensityByteArr = rs.getBytes("peakIntensity");

                if(!setMassKey) {
                    sm.setMassKey(rs.getInt("massKey"));
                    setMassKey = true;
                }


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

    public void getTableNames() throws SQLException {
        ResultSet rs = specDbConnection.getMetaData().getTables(null, null, null, null);
        while (rs.next()) {
            System.out.println(rs.getString("TABLE_NAME"));
        }
        System.out.println(specDbConnection.getSchema());
    }

    boolean isConnected(){
        return isConnected;
    }

    boolean connect() {
        try {
            // db parameters
            String url = "jdbc:sqlite:" + databaseFileName;
            // create a connection to the database
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
            specDbConnection = DriverManager.getConnection(url);//,config.toProperties());
            isConnected = true;

            //System.out.println("Connection to SQLite has been established.");

        } catch (SQLException e) {
            System.out.println(e.getMessage());
            isConnected = false;
            return false;
        }
        return true;
    }

    void disconnect() {
        try {
            if (specDbConnection != null) {
                specDbConnection.close();
                isConnected = false;
                //System.out.println("Connection to SQLite has been closed.");
            }
        } catch (SQLException ex) {
            System.out.println(ex.getMessage());
        }
    }

    void createNewMergedSpectraTable(){

        // SQL statement for creating a new table
        String sql = "CREATE TABLE IF NOT EXISTS mergedSpectraTable (\n"
                + "    peptideID integer PRIMARY KEY,\n"
                + "    peakMZ Binary,\n"
                + "    peakIntensity Binary,\n"
                + "    numElems integer\n"
                + ");";

        String sqlDrop = "DROP TABLE IF EXISTS mergedSpectraTable;";

        try (Statement stmt  = specDbConnection.createStatement()
        ){
            stmt.execute(sqlDrop);
            stmt.execute(sql);
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Insert a new row into the warehouses table
     *
     * @param peptideID
     * @param peakMzByteArr
     * @param peakIntensityByteArr
     * @param MassKey
     */
    public void insert(int peptideID, byte[] peakMzByteArr, byte[] peakIntensityByteArr, int MassKey) {
        String sql = "INSERT INTO mergedSpectraTable(peptideID,peakMz,peakIntensity,numElems) VALUES(?,?,?,?)";

        System.out.println("Peptide id: " + peptideID);

        try (PreparedStatement pstmt = specDbConnection.prepareStatement(sql)
        ){
            pstmt.setInt(1, peptideID);
            pstmt.setBytes(2, peakMzByteArr);
            pstmt.setBytes(3, peakIntensityByteArr);
            pstmt.setInt(4, MassKey);
            pstmt.executeUpdate();
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Insert a new row into the warehouses table
     *
     */
    public synchronized void insertBatch(int startPeptideID, LinkedList<float[]> listOfSpecArrays, LinkedList<Integer> numElems) throws SQLException {

        String SQL = "INSERT INTO mergedSpectraTable(peptideID,peakMz,peakIntensity,numElems) VALUES(?,?,?,?)";
        PreparedStatement pstmt = specDbConnection.prepareStatement(SQL);
        specDbConnection.setAutoCommit(false);

        for (int i = 0; i < listOfSpecArrays.size(); i++){
            float[] peakMzFloatArr;
            float[] peakIntensityFloatArr;
            byte[] peakMzByteArr;
            byte[] peakIntensityByteArr;
            peakMzFloatArr = new float[numElems.get(i)];
            peakIntensityFloatArr = new float[numElems.get(i)];

            int count = 0;
            for(int j = 0; j<listOfSpecArrays.get(i).length - 1; j++){
                if(listOfSpecArrays.get(i)[j] != 0){
                    peakMzFloatArr[count] = ((float)j)/1000;
                    peakIntensityFloatArr[count] = listOfSpecArrays.get(i)[j];
                    count++;
                }
            }

            peakMzByteArr = FloatArray2ByteArray(peakMzFloatArr);
            peakIntensityByteArr = FloatArray2ByteArray(peakIntensityFloatArr);
            int peptideID = startPeptideID + i;

            pstmt.setInt(1, peptideID);
            pstmt.setBytes(2, peakMzByteArr);
            pstmt.setBytes(3, peakIntensityByteArr);
            pstmt.setInt(4, numElems.get(i));

            // Add above SQL statement in the batch.
            pstmt.addBatch();
        }
        //Create an int[] to hold returned values
        int[] count = pstmt.executeBatch();

        //Explicitly commit statements to apply changes
        specDbConnection.commit();

    }

    void addMergedPeptideRow(int peptideID, float[] specArr, int numElems){
        float[] peakMzFloatArr;
        float[] peakIntensityFloatArr;
        byte[] peakMzByteArr;
        byte[] peakIntensityByteArr;

        peakMzFloatArr = new float[numElems];
        peakIntensityFloatArr = new float[numElems];

        int count = 0;
        for(int i = 0; i<specArr.length - 1; i++){
            if(specArr[i] != 0){
                peakMzFloatArr[count] = ((float)i)/1000;
                peakIntensityFloatArr[count] = specArr[i];
                count++;
            }
        }

        peakMzByteArr = FloatArray2ByteArray(peakMzFloatArr);
        peakIntensityByteArr = FloatArray2ByteArray(peakIntensityFloatArr);

        insert(peptideID, peakMzByteArr, peakIntensityByteArr, numElems);
    }

    public static byte[] FloatArray2ByteArray(float[] values){
        ByteBuffer buffer = ByteBuffer.allocate(4 * values.length);

        for (float value : values){
            buffer.putFloat(value);
        }

        return buffer.array();
    }


    public static void main(String[] args) throws SQLException {
    }
}
