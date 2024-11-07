// src/Main.java

import ExonSkipping.*;
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

        ArrayList<Long> gtfcds = new ArrayList<>();
        ArrayList<Long> eventsTime = new ArrayList<>();
        ArrayList<Long> write = new ArrayList<>();

        long startTime = System.currentTimeMillis();
        Genome genome = gtfParser.parseGenome(parser.getArgumentValue("-gtf"));
        long endTime = System.currentTimeMillis();
        gtfcds.add(endTime - startTime);

        startTime = System.currentTimeMillis();
        ArrayList<ExonSkippingEvent> events = genome.getExonSkippingEvents();
        endTime = System.currentTimeMillis();
        eventsTime.add(endTime - startTime);

        startTime = System.currentTimeMillis();
        Genome.writeToFile(events, parser.getArgumentValue("-o"));
        endTime = System.currentTimeMillis();
        write.add(endTime - startTime);

        //System.out.println(gtfcds.stream().mapToLong(x -> x).sum() / gtfcds.size());
        //System.out.println(eventsTime.stream().mapToLong(x -> x).sum() / gtfcds.size());
        //System.out.println(write.stream().mapToLong(x -> x).sum() / gtfcds.size());
    }
}