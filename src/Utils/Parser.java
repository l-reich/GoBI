package Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;

import ExonSkipping.*;

public class Parser {
    public Genome parseGenome(String filePath) {
        Genome genome = new Genome();
        //ArrayList<Gene> geneArrayList = genome.getGenes();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            Gene gene = null;
            Transcript transcript = null;
            RegionVector cds = null;
            String currentGene = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue; // Skip comments
                }
                GTFEntry entry = parseLine(line);
                if (!"CDS".equals(entry.getFeature())) {
                    continue; // Only process CDS lines
                }

                String geneId = entry.getAttributes().get("gene_id");
                if (geneId == null) {
                    continue; // Skip if gene_id is not present
                }
                if (!geneId.equals(currentGene)) {
                    gene = new Gene(geneId, entry.getSeqname(), entry.getStrand(), entry.getAttributes().get("gene_name"));
                    genome.addGene(gene);
                    currentGene = geneId;
                }
                /*String transId = entry.getAttributes().get("transcript_id");
                transcript = gene.getTranscripts().get(transId);
                if (transcript == null) {
                    transcript = new Transcript(transId, geneId);
                    gene.addTranscript(transcript);
                }*/

                String proteinId = entry.getAttributes().get("protein_id");
                cds = gene.getCodingSequences().get(proteinId);
                if (cds == null) {
                    cds = new RegionVector();
                    cds.addCDS(new CDS(entry.getStart(), entry.getEnd()));
                    gene.addCodingSequence(proteinId, cds);
                }
                else {
                    cds.addCDS(new CDS(entry.getStart(), entry.getEnd()));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return genome;
    }

    /*public Genome parseGenome(String filePath) {
        Genome genome = new Genome();
        //ArrayList<Gene> geneArrayList = genome.getGenes();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            Gene gene = null;
            Transcript transcript = null;
            RegionVector cds = null;
            String currentGene = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue; // Skip comments
                }
                GTFEntry entry = parseLine(line);
                if (!"CDS".equals(entry.getFeature())) {
                    continue; // Only process CDS lines
                }

                String geneId = entry.getAttributes().get("gene_id");
                if (geneId == null) {
                    continue; // Skip if gene_id is not present
                }
                if (!geneId.equals(currentGene)) {
                    gene = new Gene(geneId, entry.getSeqname(), entry.getStrand());
                    genome.addGene(gene);
                    currentGene = geneId;
                }
                String transId = entry.getAttributes().get("transcript_id");
                transcript = gene.getTranscripts().get(transId);
                if (transcript == null) {
                    transcript = new Transcript(transId, geneId);
                    gene.addTranscript(transcript);
                }
                transcript.addCDS(new CDS(entry.getStart(), entry.getEnd()));
                // Additional processing for the gene can be done here

                String proteinId = entry.getAttributes().get("protein_id");
                cds = gene.getCodingSequences().get(proteinId);
                if (cds == null) {
                    cds = new RegionVector();
                    cds.addCDS(new CDS(entry.getStart(), entry.getEnd()));
                    gene.addCodingSequence(proteinId, cds);
                }
                else {
                    cds.addCDS(new CDS(entry.getStart(), entry.getEnd()));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return genome;
    }*/

    private GTFEntry parseLine(String line) {
        String[] fields = line.split("\t");
        if (fields.length != 9) {
            throw new IllegalArgumentException("Invalid GTF format");
        }

        GTFEntry entry = new GTFEntry();
        entry.setSeqname(fields[0]);
        entry.setSource(fields[1]);
        entry.setFeature(fields[2]);
        entry.setStart(Integer.parseInt(fields[3]));
        entry.setEnd(Integer.parseInt(fields[4]));
        entry.setScore(fields[5]);
        entry.setStrand(fields[6]);
        entry.setFrame(fields[7]);
        entry.setAttributes(parseAttributes(fields[8]));

        return entry;
    }

    private Map<String, String> parseAttributes(String attributeField) {
        Map<String, String> attributes = new HashMap<>();
        String[] attributePairs = attributeField.split(";");
        for (String pair : attributePairs) {
            if (pair.trim().isEmpty()) continue;
            String[] keyValue = pair.trim().split(" ");
            if (keyValue.length == 2) {
                attributes.put(keyValue[0], keyValue[1].replaceAll("\"", ""));
            }
        }
        return attributes;
    }
}