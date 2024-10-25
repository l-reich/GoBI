package ExonSkipping;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Transcript {
    private String transcriptId;
    private String geneId;
    List<CDS> codingSequences = new ArrayList<>();
    List<Intron> introns = new ArrayList<>();

    public Transcript(String transcriptId, String geneId) {
        this.transcriptId = transcriptId;
        this.geneId = geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public void addCDS(CDS cds) {
        codingSequences.add(cds);
    }

    public List<CDS> getCodingSequences() {
        return codingSequences;
    }

    public void setCodingSequences(List<CDS> codingSequences) {
        this.codingSequences = codingSequences;
    }

    public List<Intron> getIntrons() {
        return introns;
    }

    public void setIntrons(List<Intron> introns) {
        this.introns = introns;
    }

    public void calculateIntrons() {
        introns.clear();
        if (codingSequences.size() < 2) {
            return; // No introns if there are less than 2 CDSs
        }

        // Sort CDSs by start position
        codingSequences.sort((cds1, cds2) -> Integer.compare(cds1.getStart(), cds2.getStart()));

        for (int i = 0; i < codingSequences.size() - 1; i++) {
            CDS currentCDS = codingSequences.get(i);
            CDS nextCDS = codingSequences.get(i + 1);
            int intronStart = currentCDS.getEnd() + 1;
            int intronEnd = nextCDS.getStart() - 1;
            if (intronStart <= intronEnd) {
                introns.add(new Intron(intronStart, intronEnd));
            }
        }
    }
}