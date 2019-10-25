package SpectrumMerger;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.sql.*;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;

public class spectrumDbUtil {

    private String databaseFileName;
    private Connection specDbConnection;
    private static boolean isConnected = false;

    public spectrumDbUtil(String dbFileName){
        this.databaseFileName = dbFileName;
    }

    spectrumDbUtil(){}

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

                //combine 2 arrays into peak list array
                for(int i=0;i<peakMzArray.length;i++){
                    Peak p = new Peak(peakMzArray[i], peakIntensityArray[i]);
                    list.addPeak(p);
                }
            }
            return list;
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
        return null;
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
            specDbConnection = DriverManager.getConnection(url);
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
        spectrumDbUtil db = new spectrumDbUtil("testLibDuplicateSpectraMerged.db");
        db.connect();
        db.getTableNames();
        db.createNewMergedSpectraTable();
        db.insert(5, new byte[10], new byte[10], 20000);
        db.disconnect();
    }
}
