package ReadSimulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Transcript {
    private String transcriptId;
    private List<Exon> exons;

    public Transcript(String transcriptId) {
        this.transcriptId = transcriptId;
        this.exons = new ArrayList<>();
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void addExon(Exon exon) {
        this.exons.add(exon);
    }

    public void setExons(List<Exon> exons) {
        this.exons = exons;
    }

    public List<Exon> getExons() {
        return exons;
    }
}
