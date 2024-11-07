// src/Main.java

import ReadSimulator.*;
import ReadSimulator.Gene;
import Utils.CmdParser;
import Utils.Parser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class RSMain {
    public static void main(String[] args) throws IOException {
        CmdParser parser = new CmdParser();
        parser.declareArgument("-gtf");
        parser.declareArgument("-od");
        parser.declareArgument("-length");
        parser.declareArgument("-frlength");
        parser.declareArgument("-SD");
        parser.declareArgument("-readcounts");
        parser.declareArgument("-mutationrate");
        parser.declareArgument("-fasta");
        parser.declareArgument("-fidx");

        parser.parse(args);

        RCReader rcReader = new RCReader(parser.getArgumentValue("-readcounts"));
        HashMap<String, HashMap<String, Integer>> geneTranscriptCounts = rcReader.getGeneTranscriptCounts();

        GTFParser gtfParser = new GTFParser(geneTranscriptCounts);
        RSGenome genome = gtfParser.parseGTFFile(parser.getArgumentValue("-gtf"));

        File fastaFile = new File(parser.getArgumentValue("-fasta"));
        File idxFile = new File(parser.getArgumentValue("-fidx"));
        GenomeSequenceExtractor gse = new GenomeSequenceExtractor(fastaFile, idxFile);

        System.out.println(genome.getGenes().get("ENSG00000241978").getStart() + " " + genome.getGenes().get("ENSG00000241978").getEnd());
        System.out.println(genome.getGenes().get("ENSG00000241978").getChromosome());

        String str = gse.getSequence("15", 20216406, 20216422, false);
        System.out.println(str);
    }
}