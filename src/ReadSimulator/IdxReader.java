package ReadSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

public class IdxReader {
    private HashMap<String, ArrayList<Long>> chromosomeData;

    public IdxReader(String filePath) throws IOException {
        chromosomeData = new HashMap<>();
        readIdxFile(filePath);
    }

    private void readIdxFile(String filePath) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] columns = line.split("\t");
                if (columns.length == 5) {
                    String chromosome = columns[0];
                    long start = Long.parseLong(columns[2]);
                    int charCount = Integer.parseInt(columns[3]);
                    int charCountWithNewline = Integer.parseInt(columns[4]);

                    ArrayList<Long> values = new ArrayList<>();
                    values.add(start);
                    values.add((long) charCount);
                    values.add((long) charCountWithNewline);

                    chromosomeData.put(chromosome, values);
                }
            }
        }
    }

    public HashMap<String, ArrayList<Long>> getChromosomeData() {
        return chromosomeData;
    }
}