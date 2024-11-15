package ReadSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class GTFParser {
    private HashMap<String, HashMap<String, Integer>> geneTranscriptCounts;
    private RSGenome RSGenome;

    public GTFParser(HashMap<String, HashMap<String, Integer>> geneTranscriptCounts) {
        this.geneTranscriptCounts = geneTranscriptCounts;
        this.RSGenome = new RSGenome();
    }

    public RSGenome parseGTFFile(String filePath) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            Gene currentGene = null;
            while ((line = br.readLine()) != null) {
                String[] columns = line.split("\t");
                if (columns.length < 9) continue;

                String featureType = columns[2];

                if ("gene".equals(featureType) || "exon".equals(featureType)) {
                    String attributes = columns[8];
                    String geneId = getAttribute(attributes, "gene_id");
                    String transcriptId = getAttribute(attributes, "transcript_id");

                    if ("gene".equals(featureType)) {
                        if (geneTranscriptCounts.containsKey(geneId) && !RSGenome.isGeneSaved(geneId)) {
                            currentGene = new Gene(geneId, Integer.parseInt(columns[3]), Integer.parseInt(columns[4]), columns[6], columns[0]);
                            RSGenome.addGene(geneId, currentGene);
                        }
                    } else if ("exon".equals(featureType) && currentGene != null && geneId.equals(currentGene.getGeneId())) {
                        if (transcriptId != null) {
                            if (geneTranscriptCounts.containsKey(geneId) && geneTranscriptCounts.get(geneId).containsKey(transcriptId)) {
                                Transcript transcript = getOrCreateTranscript(currentGene, transcriptId);
                                Exon exon = new Exon(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]));
                                transcript.addExon(exon);
                            }
                        }
                    }
                }
            }
        }
        return RSGenome;
    }

    private String getAttribute(String attributes, String key) {
        for (String attribute : attributes.split(";")) {
            String[] keyValue = attribute.trim().split(" ");
            if (keyValue.length == 2 && keyValue[0].equals(key)) {
                return keyValue[1].replaceAll("\"", "");
            }
        }
        return null;
    }

    private Transcript getOrCreateTranscript(Gene gene, String transcriptId) {
        for (Transcript transcript : gene.getTranscripts()) {
            if (transcript.getTranscriptId().equals(transcriptId)) {
                return transcript;
            }
        }
        Transcript newTranscript = new Transcript(transcriptId);
        gene.addTranscript(newTranscript);
        return newTranscript;
    }
}