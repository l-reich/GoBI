package ReadSimulator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

public class Writer {
    private String outputDirectory;

    public Writer(String outputDirectory) {
        this.outputDirectory = outputDirectory;
    }

    public void writeMappingInfo(List<Read> reads) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + "/read.mappinginfo"))) {
            writer.write("readid\tchr\tgene\ttranscript\tfw_regvec\trw_regvec\tt_fw_regvec\tt_rw_regvec\tfw_mut\trw_mut\n");
            for (Read read : reads) {
                String forwardGenomicRegionStr = read.getForwardGenomicRegion().stream()
                        .map(exon -> exon.getStart() + "-" + exon.getEnd())
                        .collect(Collectors.joining("|"));

                String reverseGenomicRegionStr = read.getReverseGenomicRegion().stream()
                        .map(exon -> exon.getStart() + "-" + exon.getEnd())
                        .collect(Collectors.joining("|"));

                String fwMutStr = read.getMutatedPositionsFw().stream()
                        .map(String::valueOf)
                        .collect(Collectors.joining(","));

                String rwMutStr = read.getMutatedPositionsRw().stream()
                        .map(String::valueOf)
                        .collect(Collectors.joining(","));

                writer.write(read.getReadId() + "\t" + read.getChrId() + "\t" + read.getGeneId() + "\t" + read.getTranscriptId() + "\t" +
                        forwardGenomicRegionStr + "\t" + reverseGenomicRegionStr + "\t" +
                        read.getForwardTranscriptRegion().getStart() + "-" + read.getForwardTranscriptRegion().getEnd() + "\t" +
                        read.getReverseTranscriptRegion().getStart() + "-" + read.getReverseTranscriptRegion().getEnd() + "\t" +
                        fwMutStr + "\t" + rwMutStr + "\n");
            }
        }
    }

    public void writeForwardFastq(List<Read> reads) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + "/fw.fastq"))) {
            for (Read read : reads) {
                writer.write("@" + read.getReadId() + "\n");
                writer.write(read.getForwardRead() + "\n");
                writer.write("+\n");
                writer.write("l".repeat(read.getForwardRead().length()) + "\n");
            }
        }
    }

    public void writeReverseFastq(List<Read> reads) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + "/rw.fastq"))) {
            for (Read read : reads) {
                writer.write("@" + read.getReadId() + "\n");
                writer.write(read.getReverseRead() + "\n");
                writer.write("+\n");
                writer.write("l".repeat(read.getReverseRead().length()) + "\n");
            }
        }
    }
}