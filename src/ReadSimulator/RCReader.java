package ReadSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class RCReader {
    private HashMap<String, HashMap<String, Integer>> geneTranscriptCounts;

    public RCReader(String filePath) throws IOException {
        geneTranscriptCounts = new HashMap<>();
        readTSVFile(filePath);
    }

    private void readTSVFile(String filePath) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean isFirstLine = true;
            while ((line = br.readLine()) != null) {
                if (isFirstLine) {
                    isFirstLine = false;
                    continue; // Skip the first line
                }
                String[] columns = line.split("\t");
                if (columns.length == 3) {
                    String gene = columns[0];
                    String transcript = columns[1];
                    int count = Integer.parseInt(columns[2]);

                    geneTranscriptCounts
                            .computeIfAbsent(gene, k -> new HashMap<>())
                            .put(transcript, count);
                }
            }
        }
    }

    public HashMap<String, HashMap<String, Integer>> getGeneTranscriptCounts() {
        return geneTranscriptCounts;
    }
}