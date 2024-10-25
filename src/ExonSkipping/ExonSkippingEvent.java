package ExonSkipping;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashSet;

public class ExonSkippingEvent {
    private String geneId;
    private String name;
    private String chromosome;
    private String strand;
    private int nprots;
    private Intron svIntron;
    private HashSet<Intron> wtIntrons;
    private HashSet<String> swProts;
    private HashSet<String> wtProts;
    private int minExons;
    private int maxExons;
    private int minBases;
    private int maxBases;

    public ExonSkippingEvent(String geneId, String name, String chromosome, String strand, int nprots, Intron svIntron, HashSet<Intron> wtIntrons, HashSet<String> swProts, HashSet<String> wtProts, int minExons, int maxExons, int minBases, int maxBases) {
        this.geneId = geneId;
        this.name = name;
        this.chromosome = chromosome;
        this.strand = strand;
        this.nprots = nprots;
        this.svIntron = svIntron;
        this.wtIntrons = wtIntrons;
        this.swProts = swProts;
        this.wtProts = wtProts;
        this.minExons = minExons;
        this.maxExons = maxExons;
        this.minBases = minBases;
        this.maxBases = maxBases;
    }

    @Override
    public String toString() {
        StringBuilder wtIntronsString = new StringBuilder();
        for (Intron intron : wtIntrons) {
            if (wtIntronsString.length() > 0) {
                wtIntronsString.append("|");
            }
            wtIntronsString.append(intron.getStart()).append(":").append(intron.getEnd());
        }
        return geneId + "\t" + name + "\t" + chromosome + "\t" + strand + "\t" + nprots + "\t" + nprots + "\t" +
                svIntron.getStart() + ":" + (svIntron.getEnd()+1) + "\t" + wtIntronsString.toString() + "\t" +
                String.join("|", swProts) + "\t" + String.join("|", wtProts) + "\t" + minExons + "\t" +
                maxExons + "\t" + minBases + "\t" + maxBases;
    }
}
