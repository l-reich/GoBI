// src/RSMain.java

import ReadSimulator.*;
import Utils.CmdParser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.*;
import java.util.stream.Collectors;


public class RSTest {
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

        genome.getGeneSequences(fastaFile, idxFile);

        RSGenome test = new RSGenome();
        test.addGene("ENSG00000160767", genome.getGenes().get("ENSG00000160767"));
        //test.addGene("ENSG00000232339", genome.getGenes().get("ENSG00000232339"));



        // Remove the percentage sign and convert to double
        String mutationRateStr = parser.getArgumentValue("-mutationrate").replace("%", "").trim();
        double mutationRate = Double.parseDouble(mutationRateStr) / 100.0;

        ArrayList<Read> reads = test.getReads(
                Double.parseDouble(parser.getArgumentValue("-frlength")),
                Double.parseDouble(parser.getArgumentValue("-SD")),
                Integer.parseInt(parser.getArgumentValue("-length")),
                mutationRate,
                geneTranscriptCounts);

        Writer writer = new Writer(parser.getArgumentValue("-od"));
        writer.writeForwardFastq(reads);
        writer.writeReverseFastq(reads);
        writer.writeMappingInfo(reads);
    }
}