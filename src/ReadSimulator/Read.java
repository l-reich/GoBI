package ReadSimulator;

import java.util.List;

public class Read {
    private String readId;
    private String chrId;
    private String geneId;
    private String transcriptId;
    private String forwardRead;
    private String reverseRead;
    private List<Exon> forwardGenomicRegion;
    private List<Exon> reverseGenomicRegion;
    private Exon forwardTranscriptRegion;
    private Exon reverseTranscriptRegion;
    private List<Integer> mutatedPositionsFw;
    private List<Integer> mutatedPositionsRw;
    private int fragmentLength;
    private Exon region;
    private Exon transcriptRegion;

    public Read(String readId, String chrId, String geneId, String transcriptId, String forwardRead, String reverseRead,
                List<Exon> forwardGenomicRegion, List<Exon> reverseGenomicRegion,
                Exon forwardTranscriptRegion, Exon reverseTranscriptRegion,
                List<Integer> mutatedPositionsFw, List<Integer> mutatedPositionsRw) {
        this.readId = readId;
        this.chrId = chrId;
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.forwardRead = forwardRead;
        this.reverseRead = reverseRead;
        this.forwardGenomicRegion = forwardGenomicRegion;
        this.reverseGenomicRegion = reverseGenomicRegion;
        this.forwardTranscriptRegion = forwardTranscriptRegion;
        this.reverseTranscriptRegion = reverseTranscriptRegion;
        this.mutatedPositionsFw = mutatedPositionsFw;
        this.mutatedPositionsRw = mutatedPositionsRw;
    }

    // Getters and setters for all attributes
    public String getReadId() {
        return readId;
    }

    public Exon getTranscriptRegion() {
        return transcriptRegion;
    }

    public void setTranscriptRegion(Exon transcriptRegion) {
        this.transcriptRegion = transcriptRegion;
    }

    public Exon getRegion() {
        return region;
    }

    public void setRegion(Exon region) {
        this.region = region;
    }

    public int getFragmentLength() {
        return fragmentLength;
    }

    public void setFragmentLength(int fragmentLength) {
        this.fragmentLength = fragmentLength;
    }

    public void setReadId(String readId) {
        this.readId = readId;
    }

    public String getChrId() {
        return chrId;
    }

    public void setChrId(String chrId) {
        this.chrId = chrId;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }

    public String getForwardRead() {
        return forwardRead;
    }

    public void setForwardRead(String forwardRead) {
        this.forwardRead = forwardRead;
    }

    public String getReverseRead() {
        return reverseRead;
    }

    public void setReverseRead(String reverseRead) {
        this.reverseRead = reverseRead;
    }

    public List<Exon> getForwardGenomicRegion() {
        return forwardGenomicRegion;
    }

    public void setForwardGenomicRegion(List<Exon> forwardGenomicRegion) {
        this.forwardGenomicRegion = forwardGenomicRegion;
    }

    public List<Exon> getReverseGenomicRegion() {
        return reverseGenomicRegion;
    }

    public void setReverseGenomicRegion(List<Exon> reverseGenomicRegion) {
        this.reverseGenomicRegion = reverseGenomicRegion;
    }

    public Exon getForwardTranscriptRegion() {
        return forwardTranscriptRegion;
    }

    public void setForwardTranscriptRegion(Exon forwardTranscriptRegion) {
        this.forwardTranscriptRegion = forwardTranscriptRegion;
    }

    public Exon getReverseTranscriptRegion() {
        return reverseTranscriptRegion;
    }

    public void setReverseTranscriptRegion(Exon reverseTranscriptRegion) {
        this.reverseTranscriptRegion = reverseTranscriptRegion;
    }

    public List<Integer> getMutatedPositionsFw() {
        return mutatedPositionsFw;
    }

    public void setMutatedPositionsFw(List<Integer> mutatedPositionsFw) {
        this.mutatedPositionsFw = mutatedPositionsFw;
    }

    public List<Integer> getMutatedPositionsRw() {
        return mutatedPositionsRw;
    }

    public void setMutatedPositionsRw(List<Integer> mutatedPositionsRw) {
        this.mutatedPositionsRw = mutatedPositionsRw;
    }
}