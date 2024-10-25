package ExonSkipping;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

public class Genome {
    public ArrayList<Gene> genes = new ArrayList<>();


    public ArrayList<Gene> getGenes() {
        return genes;
    }

    public void addGene(Gene gene) {
        genes.add(gene);
    }

    private void calculateGenomeIntrons() {
        for (Gene gene : genes) {
            for (Transcript transcript : gene.getTranscripts().values()) {
                transcript.calculateIntrons();
            }
        }
    }

    public void getExonSkippingEvents(String filepath) {
        ArrayList<ExonSkippingEvent> exonSkippingEvents = new ArrayList<>();
        for (Gene gene : genes) {
            ArrayList<ExonSkippingEvent> geneEvents = gene.getEvents();
            exonSkippingEvents.addAll(geneEvents);
        }
        writeToFile(exonSkippingEvents, filepath);
    }

    public static void writeToFile(ArrayList<ExonSkippingEvent> events, String filepath) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filepath))) {
            writer.write("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tSV_prots\tWT_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases");
            writer.newLine();

            for (ExonSkippingEvent event : events) {
                writer.write(event.toString());
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}