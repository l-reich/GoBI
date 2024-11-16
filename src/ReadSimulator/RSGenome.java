package ReadSimulator;

import java.io.File;
import java.io.IOException;
import java.util.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
//import org.apache.commons.math4.distribution.NormalDistribution;


public class RSGenome {
    private HashMap<String, Gene> genes;

    public RSGenome() {
        this.genes = new HashMap<>();
    }

    public void addGene(String geneId, Gene gene) {
        this.genes.put(geneId, gene);
    }

    public HashMap<String, Gene> getGenes() {
        return genes;
    }

    public boolean isGeneSaved(String geneId) {
        return genes.containsKey(geneId);
    }

    public void getGeneSequences(File fasta, File idx) throws IOException {
        GenomeSequenceExtractor gse = new GenomeSequenceExtractor(fasta, idx);
        for (String geneId : genes.keySet()) {
            Gene gene = genes.get(geneId);
            Boolean strand = gene.getStrand().equals("+");
            String sequence = gse.getSequence(gene.getChromosome(), gene.getStart(), gene.getEnd(), strand);
            gene.setSequence(sequence);
        }
    }

    public ArrayList<Read> getReads(double mean, double sd, int readLength, double mutationRate, HashMap<String, HashMap<String, Integer>> geneTranscriptCounts) {
        ArrayList<Read> reads = new ArrayList<>();
        Random random = new Random();
        int i = 0;

        for (Gene gene : genes.values()) {
            boolean isNegativeStrand = gene.getStrand().equals("-");
            if (isNegativeStrand) {
                for (Transcript transcript : gene.getTranscripts()) {
                    //Collections.reverse(transcript.getExons());
                    transcript.setExons(transcript.getExons().reversed());
                }
            }

            for (Transcript transcript : gene.getTranscripts()) {
                StringBuilder transcriptSequence = new StringBuilder();
                for (Exon exon : transcript.getExons()) {
                    String exonSequence = gene.getSequence().substring(exon.getStart() - gene.getStart(), exon.getEnd() - gene.getStart() + 1);
                    transcriptSequence.append(exonSequence);
                }
                String finalSequence;
                if (isNegativeStrand) {
                    finalSequence = getReverseComplement(transcriptSequence.toString());
                } else {
                    finalSequence = transcriptSequence.toString();
                }

                // Get the number of reads to generate for this transcript
                int numReads = geneTranscriptCounts.get(gene.getGeneId()).get(transcript.getTranscriptId());

                for (int n = 0; n < numReads; n++) {
                    int fragmentLength;
                    do {
                        fragmentLength = (int) Math.max(1, random.nextGaussian() * sd + mean);
                    } while (fragmentLength < readLength || fragmentLength > finalSequence.length());

                    int startPos = random.nextInt(finalSequence.length() - fragmentLength + 1);

                    String fragmentSequence = finalSequence.substring(startPos, startPos + fragmentLength);

                    // Generate forward and reverse reads
                    String forwardRead = fragmentSequence.substring(0, Math.min(readLength, fragmentSequence.length()));
                    String reverseRead = getReverseComplement(fragmentSequence.substring(Math.max(0, fragmentSequence.length() - readLength)));

                    // Simulate mutations
                    StringBuilder mutatedForwardRead = new StringBuilder(forwardRead);
                    StringBuilder mutatedReverseRead = new StringBuilder(reverseRead);
                    ArrayList<Integer> mutatedPositionsFw = new ArrayList<>();
                    ArrayList<Integer> mutatedPositionsRw = new ArrayList<>();

                    for (int j = 0; j < forwardRead.length(); j++) {
                        if (random.nextDouble() < mutationRate) {
                            char originalBase = forwardRead.charAt(j);
                            char mutatedBase = mutateBase(originalBase, random);
                            mutatedForwardRead.setCharAt(j, mutatedBase);
                            mutatedPositionsFw.add(j);
                        }
                    }

                    for (int j = 0; j < reverseRead.length(); j++) {
                        if (random.nextDouble() < mutationRate) {
                            char originalBase = reverseRead.charAt(j);
                            char mutatedBase = mutateBase(originalBase, random);
                            mutatedReverseRead.setCharAt(j, mutatedBase);
                            mutatedPositionsRw.add(j);
                        }
                    }
                    List<Exon> forwardGenomicRegion = new ArrayList<>();
                    List<Exon> reverseGenomicRegion = new ArrayList<>();
                    int newPos = startPos;
                    if (isNegativeStrand) {
                        newPos = finalSequence.length() - (startPos + fragmentLength);

                        //forwardGenomicRegion = getMinusFwRegion(transcript, startPos, fragmentLength, readLength);
                        //reverseGenomicRegion = getMinusRwRegion(transcript, startPos, fragmentLength, readLength, finalSequence);
                    }
                    /*List<Exon> forwardGenomicRegion = new ArrayList<>();
                    int remainingLength = readLength+1-1;
                    int currentPos = startPos;

                    for (Exon exon : transcript.getExons()) {
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        if (currentPos < exonLength) {
                            int exonStart = exon.getStart() + currentPos;
                            int exonEnd = Math.min(exon.getEnd(), exonStart + remainingLength - 1+1);
                            forwardGenomicRegion.add(new Exon(exonStart, exonEnd));
                            remainingLength -= (exonEnd - exonStart + 1);
                            if (remainingLength <= 0) break;
                            currentPos = 0;
                        } else {
                            //currentPos -= exonLength;
                        }
                    }

                    ArrayList<Exon> reverseGenomicRegion = new ArrayList<>();

                    remainingLength = readLength;
                    //currentPos = startPos + fragmentLength - 1;
                    int transcriptStart = transcript.getExons().get(0).getStart();
                    //int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
                    //currentPos = transcriptEnd -
                    currentPos = startPos + fragmentLength - 1 + transcriptStart;

                    for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                        Exon exon = transcript.getExons().get(j);
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        int exonGenomicStart = exon.getStart();
                        int exonGenomicEnd = exon.getEnd();
                        if (currentPos >= exonGenomicStart && currentPos <= exonGenomicEnd) {
                            int exonEnd = currentPos;
                            int exonStart = Math.max(exonGenomicStart, exonEnd - remainingLength + 1);
                            reverseGenomicRegion.add(new Exon(exonStart, exonEnd));
                            remainingLength -= (exonEnd - exonStart + 1);
                            if (remainingLength <= 0) break;
                            currentPos = exonGenomicStart - 1;
                        }
                    }
                     */

                    int transcriptStart = transcript.getExons().get(0).getStart();
                    int remainingLength = readLength + 1;
                    int currentPos = transcriptStart;

                    int count = newPos;
                    int gesLength = 0;
                    for (int j = 0; j < transcript.getExons().size(); j++) {
                        Exon exon = transcript.getExons().get(j);
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        gesLength += exonLength;
                            /*System.out.println(exon.getStart());
                            System.out.println(exon.getEnd());
                            System.out.println();*/
                        if (count > gesLength) {
                            //currentPos += exonLength;
                        } else {
                            currentPos = exon.getStart() + count - (gesLength - exonLength);

                            if (gesLength - count >= readLength) {
                                forwardGenomicRegion.add(new Exon(currentPos, currentPos + readLength));
                                break;
                            } else {
                                if (currentPos <= exon.getEnd()) {
                                    forwardGenomicRegion.add(new Exon(currentPos, exon.getEnd() + 1));
                                }
                                Exon next = transcript.getExons().get(j + 1);
                                int nextLength = next.getEnd() - next.getStart() + 1;
                                int remaining = readLength - (gesLength - count);

                                if (nextLength >= remaining) {
                                    forwardGenomicRegion.add(new Exon(next.getStart() - 1 + 1, next.getStart() + remaining));
                                    break;
                                } else {
                                    forwardGenomicRegion.add(new Exon(next.getStart() - 1 + 1, next.getStart() + nextLength));
                                    remaining -= nextLength;
                                    Exon next2 = transcript.getExons().get(j + 2);
                                    forwardGenomicRegion.add(new Exon(next2.getStart(), next2.getStart() + remaining));
                                    break;
                                }
                                //forwardGenomicRegion.add(new Exon(transcript.getExons().get(j + 1).getStart() - 1 + 1, transcript.getExons().get(j + 1).getStart() + readLength - (gesLength - count)));
                                //break;
                            }
                        }
                    }
                    /*int transcriptStart = transcript.getExons().get(0).getStart();
                    int remainingLength = readLength + 1;
                    int currentPos = transcriptStart;

                    int count = startPos;
                    int gesLength = 0;
                    for (Exon exon : transcript.getExons()){
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        gesLength += exonLength;
                        if (count >= gesLength){
                            //currentPos += exonLength;
                        } else {
                            currentPos = exon.getStart()+count-(gesLength-exonLength);
                            break;
                        }
                    }*/
                    // Calculate currentPosRev for reverse strand
                    int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
                    int remainingLengthRev = readLength + 1;
                    int currentPosRev = 0;
                    /*int revCount = finalSequence.length() - (startPos + fragmentLength - 1);
                    gesLength = 0;
                    for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                        Exon exon = transcript.getExons().get(j);
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        gesLength += exonLength;
                        if (revCount <= gesLength) {
                            currentPosRev = exon.getEnd() - (revCount - (gesLength - exonLength));
                            if (gesLength - revCount >= readLength) {
                                reverseGenomicRegion.add(new Exon(currentPosRev - readLength + 1, currentPosRev-1));
                                break;
                            } else {
                                if (currentPosRev >= exon.getStart()) {
                                    reverseGenomicRegion.add(new Exon(exon.getStart(), currentPosRev-1));
                                }
                                reverseGenomicRegion.add(new Exon(transcript.getExons().get(j - 1).getEnd() - readLength + 1, transcript.getExons().get(j - 1).getEnd()-1));
                                break;
                            }
                        }
                    }*/
                    int revCount = finalSequence.length() - (startPos + fragmentLength - 1 + 1);
                    if (isNegativeStrand) {
                        revCount = startPos;
                    }
                    gesLength = 0;
                    for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                        Exon exon = transcript.getExons().get(j);
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        gesLength += exonLength;
                        if (revCount < gesLength) {
                            currentPosRev = 1 + exon.getEnd() - (revCount - (gesLength - exonLength));
                            break;
                        }
                    }
                    // old fgenom calc
                    /*
                    Boolean added = false;
                    for (Exon exon : transcript.getExons()) {
                        if (exon.getEnd() >= currentPos) {
                            if (exon.getEnd() > currentPos + remainingLength - 1) {
                                if (added) {
                                    forwardGenomicRegion.add(new Exon(exon.getStart(), exon.getStart() + remainingLength - 1));
                                    break;
                                } else {
                                    forwardGenomicRegion.add(new Exon(currentPos, currentPos + remainingLength - 1));
                                    break;
                                }
                            } else {
                                forwardGenomicRegion.add(new Exon(currentPos, exon.getEnd()+1));
                                remainingLength -= (exon.getEnd() - currentPos + 1);
                                currentPos = exon.getEnd() + 1 - 1;
                                added = true;
                            }
                        }
                    }*/
                    //System.out.println(currentPosRev);

                    //int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
                    /*int remainingLengthRev = readLength + 1;
                    int currentPosRev = startPos + fragmentLength - 1 + transcriptStart;*/

                    //old rev genom calc

                    Boolean addedRev = false;

                    for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                        Exon exon = transcript.getExons().get(j);
                        if (addedRev) currentPosRev = exon.getEnd();
                        if (exon.getStart() <= currentPosRev) { //+3
                            if (exon.getStart() <= currentPosRev - remainingLengthRev + 1) {
                                if (addedRev) {
                                    if (remainingLengthRev > 0) {
                                        reverseGenomicRegion.add(new Exon(exon.getEnd() + 1 + 1 - remainingLengthRev, exon.getEnd() + 1));
                                        //System.out.println(1 + " " + (exon.getEnd() + 2 - remainingLengthRev) + " " + (exon.getEnd() + 1));
                                        break;
                                    }
                                } else {
                                    // reverseGenomicRegion.add(new Exon(currentPosRev - remainingLengthRev, currentPosRev));
                                    reverseGenomicRegion.add(new Exon(currentPosRev + 1 - remainingLengthRev, currentPosRev));
                                    //System.out.println(2 + " " + (currentPosRev + 3 - remainingLengthRev) + " " + (currentPosRev + 2));
                                    break;
                                }
                            } else {
                                reverseGenomicRegion.add(new Exon(exon.getStart(), currentPosRev + 1 - 1));
                                remainingLengthRev -= (currentPosRev - exon.getStart());
                                //remainingLengthRev -= (currentPosRev - exon.getEnd());
                                //currentPosRev = exon.getStart() - 1;
                                //currentPosRev = transcript.getExons().get(j - 1).getEnd();
                                addedRev = true;
                                // System.out.println(3 + " " + exon.getStart() + " " +( currentPosRev +1));
                                //System.out.println();
                            }
                        }
                    }
                    if (!isNegativeStrand) {
                        reverseGenomicRegion = reverseGenomicRegion.reversed();
                        if (reverseGenomicRegion.size() == 3) {
                            Exon first = reverseGenomicRegion.get(0);
                            Exon second = reverseGenomicRegion.get(1);
                            first.setStart(first.getStart() + 1);
                            second.setEnd(second.getEnd() + 1);
                            if (first.getStart() == first.getEnd()) {
                                reverseGenomicRegion.remove(0);
                            }
                        }
                    }

                    if (isNegativeStrand) {
                        List<Exon> temp = forwardGenomicRegion;
                        forwardGenomicRegion = reverseGenomicRegion.reversed();
                        reverseGenomicRegion = temp;
                        if (forwardGenomicRegion.size() == 3) {
                            Exon first = forwardGenomicRegion.get(0);
                            Exon second = forwardGenomicRegion.get(1);
                            first.setStart(first.getStart() + 1);
                            second.setEnd(second.getEnd() + 1);
                            if (first.getStart() == first.getEnd()) {
                                forwardGenomicRegion.remove(0);
                            }
                        }
                    }

                    Exon forwardTranscriptRegion = new Exon(startPos, startPos + readLength);
                    Exon reverseTranscriptRegion = new Exon(startPos + fragmentLength - readLength, startPos + fragmentLength);

                    Read read = new Read(
                            String.valueOf(i), // Replace with actual read ID
                            gene.getChromosome(),
                            gene.getGeneId(),
                            transcript.getTranscriptId(),
                            mutatedForwardRead.toString(), // forwardRead
                            mutatedReverseRead.toString(), // reverseRead
                            forwardGenomicRegion, // forwardGenomicRegion
                            reverseGenomicRegion, // reverseGenomicRegion
                            forwardTranscriptRegion,
                            reverseTranscriptRegion,
                            mutatedPositionsFw, // mutatedPositionsFw
                            mutatedPositionsRw // mutatedPositionsRw
                    );
                    /*read.setFragmentLength(fragmentLength);
                    read.setRegion(new Exon(startPos + transcriptStart, startPos + +transcriptStart + fragmentLength - 1));
                    read.setTranscriptRegion(new Exon(startPos, startPos + readLength - 1));*/
                    reads.add(read);
                    i++;
                }
            }
        }
        return reads;
    }

    private char mutateBase(char base, Random random) {
        char[] bases = {'A', 'T', 'C', 'G'};
        char newBase;
        do {
            newBase = bases[random.nextInt(bases.length)];
        } while (newBase == base);
        return newBase;
    }

    public List<Exon> getMinusFwRegion(Transcript transcript, int startPos, int fragmentLength, int readLength) {
        transcript.setExons(transcript.getExons().reversed());

        List<Exon> forwardGenomicRegion = new ArrayList<>();
        int remainingLength = readLength;
        int currentPos;

        for (Exon exon : transcript.getExons()) {
            System.out.println(exon.getStart());
            System.out.println(exon.getEnd());
            System.out.println();
        }

        int transcriptLength = 0;
        for (Exon exon : transcript.getExons()) {
            transcriptLength += exon.getEnd() - exon.getStart() + 1;
        }
        int count = startPos;
        int gesLength = 0;

        for (int j = 0; j < transcript.getExons().size(); j++) {
            Exon exon = transcript.getExons().get(j);
            int exonLength = exon.getEnd() - exon.getStart() + 1;
            gesLength += exonLength;

            if (count > gesLength) {
                //currentPos += exonLength;
            } else {
                currentPos = exon.getEnd() - (count - (gesLength - exonLength));
                //currentPos = exon.getStart() + (count - (gesLength - exonLength));
                System.out.println(exon.getStart());
                System.out.println(currentPos);
                System.out.println(exon.getEnd());
                System.out.println();
                System.out.println("------");

                if (gesLength - count >= readLength) {
                    forwardGenomicRegion.add(new Exon(currentPos, currentPos + readLength));
                    break;
                } else {
                   /*if (currentPos <= exon.getEnd()) { //getEnd
                        forwardGenomicRegion.add(new Exon(currentPos, currentPos + currentPos - exon.getStart()));
                        remainingLength -= (currentPos - exon.getStart());
                    }*/
                    forwardGenomicRegion.add(new Exon(currentPos, currentPos + currentPos - exon.getStart()));
                    remainingLength -= (1 + currentPos - exon.getStart());
                    Exon next = transcript.getExons().get(j + 1);
                    int gap = exon.getStart() - next.getEnd();
                    currentPos += currentPos + currentPos - exon.getStart() + gap;
                    if ((next.getEnd() - next.getStart() + 1) >= remainingLength) {
                        forwardGenomicRegion.add(new Exon(currentPos, currentPos + remainingLength));
                        //System.out.println(2 + " " + (next.getStart() - 1 + 1) + " " + (next.getStart() + readLength - (gesLength - count)));
                        break;
                    } else {
                        forwardGenomicRegion.add(new Exon(currentPos, currentPos + next.getEnd() - next.getStart() + 1));
                        remainingLength -= (next.getEnd() - next.getStart() + 1);
                        Exon next2 = transcript.getExons().get(j + 2);
                        forwardGenomicRegion.add(new Exon(next2.getStart(), next2.getStart() + remainingLength));
                        break;
                    }
                }
            }
        }

        /*for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
            Exon exon = transcript.getExons().get(j);
            int exonLength = exon.getEnd() - exon.getStart() + 1;
            if (currentPos >= exon.getStart() && currentPos <= exon.getEnd()) {
                int exonEnd = currentPos;
                int exonStart = Math.max(exon.getStart(), exonEnd - remainingLength + 1);
                forwardGenomicRegion.add(new Exon(exonStart, exonEnd + 1));
                remainingLength -= (exonEnd - exonStart + 1);
                if (remainingLength <= 0) break;
                currentPos = exon.getStart() - 1;
            }
        }*/
        return forwardGenomicRegion;
    }

    public List<Exon> getMinusFwRegionNew(Transcript transcript, int startPos, int fragmentLength, int readLength) {
        transcript.setExons(transcript.getExons().reversed());

        List<Exon> forwardGenomicRegion = new ArrayList<>();
        int remainingLength = readLength;
        int currentPos = startPos + fragmentLength - 1;

        int transcriptStart = transcript.getExons().get(0).getEnd();
        //int currentPos = transcriptStart;

        int count = startPos;
        int gesLength = 0;
        for (int j = 0; j < transcript.getExons().size(); j++) {
            Exon exon = transcript.getExons().get(j);
            int exonLength = exon.getEnd() - exon.getStart() + 1;
            gesLength += exonLength;
            if (count > gesLength) {
                //currentPos += exonLength;
            } else {
                currentPos = exon.getEnd() - (count - (gesLength - exonLength));

                if (gesLength - count >= readLength) {
                    forwardGenomicRegion.add(new Exon(currentPos, currentPos + readLength));
                    break;
                } else {
                    if (currentPos <= exon.getEnd()) {
                        forwardGenomicRegion.add(new Exon(currentPos, exon.getEnd() + 1));
                    }
                    forwardGenomicRegion.add(new Exon(transcript.getExons().get(j + 1).getStart() - 1 + 1, transcript.getExons().get(j + 1).getStart() + readLength - (gesLength - count)));
                    break;
                }
            }
        }

        /*for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
            Exon exon = transcript.getExons().get(j);
            int exonLength = exon.getEnd() - exon.getStart() + 1;
            if (currentPos >= exon.getStart() && currentPos <= exon.getEnd()) {
                int exonEnd = currentPos;
                int exonStart = Math.max(exon.getStart(), exonEnd - remainingLength + 1);
                forwardGenomicRegion.add(new Exon(exonStart, exonEnd + 1));
                remainingLength -= (exonEnd - exonStart + 1);
                if (remainingLength <= 0) break;
                currentPos = exon.getStart() - 1;
            }
        }*/
        return forwardGenomicRegion;
    }

    public List<Exon> getMinusRwRegion(Transcript transcript, int startPos, int fragmentLength, int readLength, String finalSequence) {
        List<Exon> reverseGenomicRegion = new ArrayList<>();
        // Calculate currentPosRev for reverse strand
        int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
        int remainingLengthRev = readLength + 1;
        int currentPosRev = 0;


                    /*int revCount = finalSequence.length() - (startPos + fragmentLength - 1);
                    gesLength = 0;
                    for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                        Exon exon = transcript.getExons().get(j);
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        gesLength += exonLength;
                        if (revCount <= gesLength) {
                            currentPosRev = exon.getEnd() - (revCount - (gesLength - exonLength));
                            if (gesLength - revCount >= readLength) {
                                reverseGenomicRegion.add(new Exon(currentPosRev - readLength + 1, currentPosRev-1));
                                break;
                            } else {
                                if (currentPosRev >= exon.getStart()) {
                                    reverseGenomicRegion.add(new Exon(exon.getStart(), currentPosRev-1));
                                }
                                reverseGenomicRegion.add(new Exon(transcript.getExons().get(j - 1).getEnd() - readLength + 1, transcript.getExons().get(j - 1).getEnd()-1));
                                break;
                            }
                        }
                    }*/
        int revCount = finalSequence.length() - (startPos + fragmentLength - 1 + 1);
        int gesLength = 0;
        for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
            Exon exon = transcript.getExons().get(j);
            int exonLength = exon.getEnd() - exon.getStart() + 1;
            gesLength += exonLength;
            if (revCount < gesLength) {
                currentPosRev = 1 + exon.getEnd() - (revCount - (gesLength - exonLength));
                break;
            }
        }
        // old fgenom calc
                    /*
                    Boolean added = false;
                    for (Exon exon : transcript.getExons()) {
                        if (exon.getEnd() >= currentPos) {
                            if (exon.getEnd() > currentPos + remainingLength - 1) {
                                if (added) {
                                    forwardGenomicRegion.add(new Exon(exon.getStart(), exon.getStart() + remainingLength - 1));
                                    break;
                                } else {
                                    forwardGenomicRegion.add(new Exon(currentPos, currentPos + remainingLength - 1));
                                    break;
                                }
                            } else {
                                forwardGenomicRegion.add(new Exon(currentPos, exon.getEnd()+1));
                                remainingLength -= (exon.getEnd() - currentPos + 1);
                                currentPos = exon.getEnd() + 1 - 1;
                                added = true;
                            }
                        }
                    }*/
        //System.out.println(currentPosRev);

        //int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
                    /*int remainingLengthRev = readLength + 1;
                    int currentPosRev = startPos + fragmentLength - 1 + transcriptStart;*/

        //old rev genom calc

        Boolean addedRev = false;

        for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
            Exon exon = transcript.getExons().get(j);
            if (addedRev) currentPosRev = exon.getEnd();
            if (exon.getStart() <= currentPosRev) { //+3
                if (exon.getStart() <= currentPosRev - remainingLengthRev + 1) {
                    if (addedRev) {
                        if (remainingLengthRev > 0) {
                            reverseGenomicRegion.add(new Exon(exon.getEnd() + 1 + 1 - remainingLengthRev, exon.getEnd() + 1));
                            //System.out.println(1 + " " + (exon.getEnd() + 2 - remainingLengthRev) + " " + (exon.getEnd() + 1));
                            break;
                        }
                    } else {
                        // reverseGenomicRegion.add(new Exon(currentPosRev - remainingLengthRev, currentPosRev));
                        reverseGenomicRegion.add(new Exon(currentPosRev + 1 - remainingLengthRev, currentPosRev));
                        //System.out.println(2 + " " + (currentPosRev + 3 - remainingLengthRev) + " " + (currentPosRev + 2));
                        break;
                    }
                } else {
                    reverseGenomicRegion.add(new Exon(exon.getStart(), currentPosRev + 1 - 1));
                    remainingLengthRev -= (currentPosRev - exon.getStart());
                    //remainingLengthRev -= (currentPosRev - exon.getEnd());
                    //currentPosRev = exon.getStart() - 1;
                    //currentPosRev = transcript.getExons().get(j - 1).getEnd();
                    addedRev = true;
                    // System.out.println(3 + " " + exon.getStart() + " " +( currentPosRev +1));
                    //System.out.println();
                }
            }
        }
        reverseGenomicRegion = reverseGenomicRegion.reversed();
        return reverseGenomicRegion;
    }

    private String getReverseComplement(String sequence) {
        StringBuilder reverseComplement = new StringBuilder();
        for (int i = sequence.length() - 1; i >= 0; i--) {
            char base = sequence.charAt(i);
            switch (base) {
                case 'A':
                    reverseComplement.append('T');
                    break;
                case 'T':
                    reverseComplement.append('A');
                    break;
                case 'C':
                    reverseComplement.append('G');
                    break;
                case 'G':
                    reverseComplement.append('C');
                    break;
                default:
                    reverseComplement.append(base);
                    break;
            }
        }
        return reverseComplement.toString();
    }

    private List<Exon> mapToPlus(Transcript transcript) {
        List<Exon> plusExons = new ArrayList<>();

        int transcriptStart = transcript.getExons().get(0).getStart();
        int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
        int transcriptLength = transcriptEnd - transcriptStart + 1;

        for (Exon exon : transcript.getExons()) {
            int exonStart = exon.getStart();
            int exonEnd = exon.getEnd();
            int newEnd = transcriptEnd - (exonStart - transcriptStart);
            int newStart = newEnd - (exonEnd - exonStart);
            plusExons.add(new Exon(newStart, newEnd));
        }
        return plusExons.reversed();
    }

    public List<Read> getReadsForGene(int id, Gene gene, double mean, double sd, int readLength, double mutationRate, HashMap<String, HashMap<String, Integer>> geneTranscriptCounts) {

        ArrayList<Read> reads = new ArrayList<>();
        Random random = new Random();
        int i = id;

        boolean isNegativeStrand = gene.getStrand().equals("-");
        if (isNegativeStrand) {
            for (Transcript transcript : gene.getTranscripts()) {
                //Collections.reverse(transcript.getExons());
                transcript.setExons(transcript.getExons().reversed());
            }
        }

        for (Transcript transcript : gene.getTranscripts()) {
            StringBuilder transcriptSequence = new StringBuilder();
            for (Exon exon : transcript.getExons()) {
                String exonSequence = gene.getSequence().substring(exon.getStart() - gene.getStart(), exon.getEnd() - gene.getStart() + 1);
                transcriptSequence.append(exonSequence);
            }
            String finalSequence;
            if (isNegativeStrand) {
                finalSequence = getReverseComplement(transcriptSequence.toString());
            } else {
                finalSequence = transcriptSequence.toString();
            }

            // Get the number of reads to generate for this transcript
            int numReads = geneTranscriptCounts.get(gene.getGeneId()).get(transcript.getTranscriptId());

            for (int n = 0; n < numReads; n++) {
                int fragmentLength;
                do {
                    fragmentLength = (int) Math.max(1, random.nextGaussian() * sd + mean);
                } while (fragmentLength < readLength || fragmentLength > finalSequence.length());

                int startPos = random.nextInt(finalSequence.length() - fragmentLength + 1);

                String fragmentSequence = finalSequence.substring(startPos, startPos + fragmentLength);

                // Generate forward and reverse reads
                String forwardRead = fragmentSequence.substring(0, Math.min(readLength, fragmentSequence.length()));
                String reverseRead = getReverseComplement(fragmentSequence.substring(Math.max(0, fragmentSequence.length() - readLength)));

                // Simulate mutations
                StringBuilder mutatedForwardRead = new StringBuilder(forwardRead);
                StringBuilder mutatedReverseRead = new StringBuilder(reverseRead);
                ArrayList<Integer> mutatedPositionsFw = new ArrayList<>();
                ArrayList<Integer> mutatedPositionsRw = new ArrayList<>();

                for (int j = 0; j < forwardRead.length(); j++) {
                    if (random.nextDouble() < mutationRate) {
                        char originalBase = forwardRead.charAt(j);
                        char mutatedBase = mutateBase(originalBase, random);
                        mutatedForwardRead.setCharAt(j, mutatedBase);
                        mutatedPositionsFw.add(j);
                    }
                }

                for (int j = 0; j < reverseRead.length(); j++) {
                    if (random.nextDouble() < mutationRate) {
                        char originalBase = reverseRead.charAt(j);
                        char mutatedBase = mutateBase(originalBase, random);
                        mutatedReverseRead.setCharAt(j, mutatedBase);
                        mutatedPositionsRw.add(j);
                    }
                }
                List<Exon> forwardGenomicRegion = new ArrayList<>();
                List<Exon> reverseGenomicRegion = new ArrayList<>();
                int newPos = startPos;

                if (isNegativeStrand) {
                    newPos = finalSequence.length() - (startPos + fragmentLength);
                }

                int transcriptStart = transcript.getExons().get(0).getStart();
                int remainingLength = readLength + 1;
                int currentPos = transcriptStart;

                int count = newPos;
                int gesLength = 0;
                for (int j = 0; j < transcript.getExons().size(); j++) {
                    Exon exon = transcript.getExons().get(j);
                    int exonLength = exon.getEnd() - exon.getStart() + 1;
                    gesLength += exonLength;

                    if (count > gesLength) {
                    } else {
                        currentPos = exon.getStart() + count - (gesLength - exonLength);

                        if (gesLength - count >= readLength) {
                            forwardGenomicRegion.add(new Exon(currentPos, currentPos + readLength));
                            break;
                        } else {
                            if (currentPos <= exon.getEnd()) {
                                forwardGenomicRegion.add(new Exon(currentPos, exon.getEnd() + 1));
                            }
                            Exon next = transcript.getExons().get(j + 1);
                            int nextLength = next.getEnd() - next.getStart() + 1;
                            int remaining = readLength - (gesLength - count);

                            if (nextLength >= remaining) {
                                forwardGenomicRegion.add(new Exon(next.getStart() - 1 + 1, next.getStart() + remaining));
                                break;
                            } else {
                                forwardGenomicRegion.add(new Exon(next.getStart() - 1 + 1, next.getStart() + nextLength));
                                remaining -= nextLength;
                                Exon next2 = transcript.getExons().get(j + 2);
                                forwardGenomicRegion.add(new Exon(next2.getStart(), next2.getStart() + remaining));
                                break;
                            }
                        }
                    }
                }

                int transcriptEnd = transcript.getExons().get(transcript.getExons().size() - 1).getEnd();
                int remainingLengthRev = readLength + 1;
                int currentPosRev = 0;
                    /*int revCount = finalSequence.length() - (startPos + fragmentLength - 1);
                    gesLength = 0;
                    for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                        Exon exon = transcript.getExons().get(j);
                        int exonLength = exon.getEnd() - exon.getStart() + 1;
                        gesLength += exonLength;
                        if (revCount <= gesLength) {
                            currentPosRev = exon.getEnd() - (revCount - (gesLength - exonLength));
                            if (gesLength - revCount >= readLength) {
                                reverseGenomicRegion.add(new Exon(currentPosRev - readLength + 1, currentPosRev-1));
                                break;
                            } else {
                                if (currentPosRev >= exon.getStart()) {
                                    reverseGenomicRegion.add(new Exon(exon.getStart(), currentPosRev-1));
                                }
                                reverseGenomicRegion.add(new Exon(transcript.getExons().get(j - 1).getEnd() - readLength + 1, transcript.getExons().get(j - 1).getEnd()-1));
                                break;
                            }
                        }
                    }*/
                int revCount = finalSequence.length() - (startPos + fragmentLength - 1 + 1);
                if (isNegativeStrand) {
                    revCount = startPos;
                }
                gesLength = 0;
                for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                    Exon exon = transcript.getExons().get(j);
                    int exonLength = exon.getEnd() - exon.getStart() + 1;
                    gesLength += exonLength;
                    if (revCount < gesLength) {
                        currentPosRev = 1 + exon.getEnd() - (revCount - (gesLength - exonLength));
                        break;
                    }
                }

                Boolean addedRev = false;

                for (int j = transcript.getExons().size() - 1; j >= 0; j--) {
                    Exon exon = transcript.getExons().get(j);
                    if (addedRev) currentPosRev = exon.getEnd();
                    if (exon.getStart() <= currentPosRev) { //+3
                        if (exon.getStart() <= currentPosRev - remainingLengthRev + 1) {
                            if (addedRev) {
                                if (remainingLengthRev > 0) {
                                    reverseGenomicRegion.add(new Exon(exon.getEnd() + 1 + 1 - remainingLengthRev, exon.getEnd() + 1));
                                    break;
                                }
                            } else {
                                reverseGenomicRegion.add(new Exon(currentPosRev + 1 - remainingLengthRev, currentPosRev));
                                break;
                            }
                        } else {
                            reverseGenomicRegion.add(new Exon(exon.getStart(), currentPosRev + 1 - 1));
                            remainingLengthRev -= (currentPosRev - exon.getStart());
                            addedRev = true;
                        }
                    }
                }
                if (!isNegativeStrand) {
                    reverseGenomicRegion = reverseGenomicRegion.reversed();
                    if (reverseGenomicRegion.size() == 3) {
                        Exon first = reverseGenomicRegion.get(0);
                        Exon second = reverseGenomicRegion.get(1);
                        first.setStart(first.getStart() + 1);
                        second.setEnd(second.getEnd() + 1);
                        if (first.getStart() == first.getEnd()) {
                            reverseGenomicRegion.remove(0);
                        }
                    }
                }

                if (isNegativeStrand) {
                    List<Exon> temp = forwardGenomicRegion;
                    forwardGenomicRegion = reverseGenomicRegion.reversed();
                    reverseGenomicRegion = temp;
                    if (forwardGenomicRegion.size() == 3) {
                        Exon first = forwardGenomicRegion.get(0);
                        Exon second = forwardGenomicRegion.get(1);
                        first.setStart(first.getStart() + 1);
                        second.setEnd(second.getEnd() + 1);
                        if (first.getStart() == first.getEnd()) {
                            forwardGenomicRegion.remove(0);
                        }
                    }
                }

                Exon forwardTranscriptRegion = new Exon(startPos, startPos + readLength);
                Exon reverseTranscriptRegion = new Exon(startPos + fragmentLength - readLength, startPos + fragmentLength);

                Read read = new Read(
                        String.valueOf(i), // Replace with actual read ID
                        gene.getChromosome(),
                        gene.getGeneId(),
                        transcript.getTranscriptId(),
                        mutatedForwardRead.toString(), // forwardRead
                        mutatedReverseRead.toString(), // reverseRead
                        forwardGenomicRegion, // forwardGenomicRegion
                        reverseGenomicRegion, // reverseGenomicRegion
                        forwardTranscriptRegion,
                        reverseTranscriptRegion,
                        mutatedPositionsFw, // mutatedPositionsFw
                        mutatedPositionsRw // mutatedPositionsRw
                );
                reads.add(read);
                i++;
            }
        }
        return reads;
    }
}