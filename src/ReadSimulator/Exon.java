package ReadSimulator;

public class Exon {
    private int start;
    private int end;

    public Exon(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }
}