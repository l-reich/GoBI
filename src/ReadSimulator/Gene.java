package ReadSimulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
public class Gene {
    private String geneId;
    private int start;
    private int end;
    private String strand;
    private String chromosome;
    private String sequence;
    private List<ReadSimulator.Transcript> transcripts;

    public Gene(String geneId, int start, int end, String strand, String chromosome) {
        this.geneId = geneId;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.chromosome = chromosome;
        this.transcripts = new ArrayList<>();
    }

    public String getGeneId() {
        return geneId;
    }

    public void addTranscript(ReadSimulator.Transcript transcript) {
        this.transcripts.add(transcript);
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public List<ReadSimulator.Transcript> getTranscripts() {
        return transcripts;
    }
}