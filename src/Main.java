// src/Main.java

import ExonSkipping.CDS;
import ExonSkipping.Gene;
import ExonSkipping.Genome;
import ExonSkipping.Transcript;
import Utils.CmdParser;
import Utils.Parser;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Main {
    public static void main(String[] args) {
        CmdParser parser = new CmdParser();
        parser.declareArgument("-gtf");
        parser.declareArgument("-o");

        parser.parse(args);
        Parser gtfParser = new Parser();
        Genome genome = gtfParser.parseGenome(parser.getArgumentValue("-gtf"));
        genome.getExonSkippingEvents(parser.getArgumentValue("-o"));
        /*for (Map.Entry<String, Transcript> entry : transcripts.entrySet()) {
            Transcript transcript = entry.getValue();
            //System.out.println("Transcript ID: " + transcript.getTranscriptId());
            for (int i = 0; i < transcript.getCodingSequences().size(); i++) {
                System.out.println("CDS " + i + ": " + transcript.getCodingSequences().get(i).getStart() + " - " + transcript.getCodingSequences().get(i).getEnd());
            }
        }*/

    }
}