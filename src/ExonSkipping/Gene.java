// src/ExonSkipping/Gene.java
package ExonSkipping;

import java.util.*;

public class Gene {
    private String geneId;
    private String chromosome;
    private String strand;
    private String geneName;
    private HashMap<String, Transcript> transcripts = new HashMap<>();
    private HashMap<String, RegionVector> codingSequences = new HashMap<>();
    private HashSet<Intron> introns = new HashSet<>();

    public Gene(String geneId, String chromosome, String strand, String geneName) {
        this.geneId = geneId;
        this.chromosome = chromosome;
        this.strand = strand;
        this.geneName = geneName;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
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

    public HashMap<String, Transcript> getTranscripts() {
        return transcripts;
    }

    public void setTranscripts(HashMap<String, Transcript> transcripts) {
        this.transcripts = transcripts;
    }

    public void addTranscript(Transcript transcript) {
        transcripts.put(transcript.getTranscriptId(), transcript);
    }

    public HashMap<String, RegionVector> getCodingSequences() {
        return codingSequences;
    }

    public void addCodingSequence(String proteinId, RegionVector cds) {
        codingSequences.put(proteinId, cds);
    }

    public void calculateIntrons() {
        introns.clear();
        for (RegionVector regionVector : codingSequences.values()) {
            List<CDS> codingSequences = regionVector.getCdsList();
            if (codingSequences.size() < 2) {
                continue; // No introns if there are less than 2 CDSs
            }

            // Sort CDSs by start position
            //codingSequences.sort(Comparator.comparingInt(CDS::getStart));

            for (int i = 0; i < codingSequences.size() - 1; i++) {
                CDS currentCDS = codingSequences.get(i);
                CDS nextCDS = codingSequences.get(i + 1);
                int intronStart = currentCDS.getEnd() + 1;
                int intronEnd = nextCDS.getStart() - 1;
                if (intronStart <= intronEnd) {
                    introns.add(new Intron(intronStart, intronEnd));
                    //System.out.println(geneId + "\t" + intronStart + "\t" + intronEnd);
                }
            }
        }
    }

    public ArrayList<ExonSkippingEvent> getEvents() {
        calculateIntrons();
        ArrayList<ExonSkippingEvent> exonSkippingEvents = new ArrayList<>();
        for (Intron intron : introns) {
            HashSet<String> sv = calculateSV(intron);
            HashSet<String> wtStart = calculatewtStart(intron);
            HashSet<String> wtEnd = calculatewtEnd(intron);
            wtStart.retainAll(wtEnd); // Intersection of wtStart and wtEnd
            wtStart.removeAll(sv); // Remove elements in sv

            if (!(wtStart.isEmpty())) {
                HashSet<Intron> wtIntrons = new HashSet<>();
                //System.out.println(geneId + "\t" + intron.getStart() + "\t" + intron.getEnd() + "\t" + sv.size() + "\t" + wtStart.size());
                int maxE = 0;
                int minE = Integer.MAX_VALUE;
                int maxB = 0;
                int minB = Integer.MAX_VALUE;
                for (String id : wtStart) {
                    int e = 0;
                    int b = 0;
                    boolean count = false;
                    List<CDS> cdsList = codingSequences.get(id).getCdsList();
                    for (int i = 0; i < cdsList.size() - 1; i++) {
                        CDS currentCDS = cdsList.get(i);
                        CDS nextCDS = cdsList.get(i + 1);
                        if (currentCDS.getEnd() + 1 == intron.getStart()) {
                            wtIntrons.add(new Intron(currentCDS.getEnd() + 1, nextCDS.getStart()));
                            count = true;
                            continue;
                        }
                        if (count) {
                            if (currentCDS.getStart() - 1 != intron.getEnd()) {
                                e++;
                                b += (currentCDS.getEnd() - currentCDS.getStart() + 1);
                                wtIntrons.add(new Intron(currentCDS.getEnd() + 1, nextCDS.getStart()));
                            } else {
                                break;
                            }
                        }
                    }
                    if (e > maxE) {
                        maxE = e;
                    }
                    if (e < minE) {
                        minE = e;
                    }
                    if (b > maxB) {
                        maxB = b;
                    }
                    if (b < minB) {
                        minB = b;
                    }
                }
                exonSkippingEvents.add(new ExonSkippingEvent(geneId, geneName, chromosome, strand, codingSequences.size(), intron, wtIntrons, sv, wtStart, maxE, minE, maxB, minB));
            }

        }
        return exonSkippingEvents;
    }

    /*public ArrayList<ExonSkippingEvent> getEvents() {
        calculateIntrons();
        ArrayList<ExonSkippingEvent> exonSkippingEvents = new ArrayList<>();
        for (Intron intron : introns) {
            HashSet<String> sv = calculateSV(intron);
            HashSet<String> wtStart = calculatewtStart(intron);
            HashSet<String> wtEnd = calculatewtEnd(intron);
            wtStart.retainAll(wtEnd); // Intersection of wtStart and wtEnd
            wtStart.removeAll(sv); // Remove elements in sv

            if (!(wtStart.isEmpty())) {
                HashSet<Intron> wtIntrons = new HashSet<>();
                //System.out.println(geneId + "\t" + intron.getStart() + "\t" + intron.getEnd() + "\t" + sv.size() + "\t" + wtStart.size());
                int maxE = 0; int minE = Integer.MAX_VALUE; int maxB = 0; int minB = Integer.MAX_VALUE;
                for (String id : wtStart) {
                    int e = 0; int b = 0;
                    boolean count = false;
                    List<CDS> cdsList = codingSequences.get(id).getCdsList();
                    for (CDS cds : cdsList) {
                        if (cds.getEnd()+1 == intron.getStart()) {
                            count = true;
                            continue;
                        }
                        if (count){
                            if (cds.getStart()-1 != intron.getEnd()) {

                                e++;
                                b = b + (cds.getEnd() - cds.getStart() + 1);
                            }
                            else {
                                break;
                            }
                        }
                    }
                    if (e > maxE) {
                        maxE = e;
                    }
                    if (e < minE) {
                        minE = e;
                    }
                    if (b > maxB) {
                        maxB = b;
                    }
                    if (b < minB) {
                        minB = b;
                    }
                }
                System.out.println(maxE+"\t"+minE+"\t"+maxB+"\t"+minB);
                //exonSkippingEvents.add(new ExonSkippingEvent(geneId, geneName, chromosome, strand, codingSequences.size(), intron, wtStart, maxE, minE, maxB, minB));
            }

        }


        return null;
    }*/

    private HashSet<String> calculateSV(Intron intron) {
        HashSet<String> sv = new HashSet<>();
        for (String pro_id : codingSequences.keySet()) {
            List<CDS> cdsList = codingSequences.get(pro_id).getCdsList();
            for (int i = 0; i < cdsList.size(); i++) {
                CDS cds = cdsList.get(i);
                if (cds.getEnd() + 1 == intron.getStart() && cdsList.size() > i + 1) {
                    if (cdsList.get(i + 1).getStart() - 1 == intron.getEnd()) {
                        sv.add(pro_id);
                    }
                }
            }
        }
        return sv;
    }

    private HashSet<String> calculatewtStart(Intron intron) {
        HashSet<String> wtStart = new HashSet<>();
        for (String pro_id : codingSequences.keySet()) {
            List<CDS> cdsList = codingSequences.get(pro_id).getCdsList();
            for (CDS cds : cdsList) {
                if (intron.getStart() == cds.getEnd() + 1) {
                    wtStart.add(pro_id);
                }
            }
        }
        return wtStart;
    }

    private HashSet<String> calculatewtEnd(Intron intron) {
        HashSet<String> wtEnd = new HashSet<>();
        for (String pro_id : codingSequences.keySet()) {
            List<CDS> cdsList = codingSequences.get(pro_id).getCdsList();
            for (CDS cds : cdsList) {
                if (intron.getEnd() == cds.getStart() - 1) {
                    wtEnd.add(pro_id);
                }
            }
        }
        return wtEnd;
    }

    public void setCodingSequences(HashMap<String, RegionVector> codingSequences) {
        this.codingSequences = codingSequences;
    }

    public HashSet<Intron> getIntrons() {
        return introns;
    }

    public void setIntrons(HashSet<Intron> introns) {
        this.introns = introns;
    }
}