package ReadSimulator;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;

public class GenomeSequenceExtractor {
    private RandomAccessFile raf;
    private HashMap<String, ArrayList<Long>> chromosomeData;

    public GenomeSequenceExtractor(File fasta, File idx) throws IOException {
        this.raf = new RandomAccessFile(fasta, "r");
        IdxReader idxReader = new IdxReader(idx.getPath());
        this.chromosomeData = idxReader.getChromosomeData();
    }

    public String getSequence(String chr, int start, int end, Boolean strand) throws IOException {
        if (!chromosomeData.containsKey(chr)) {
            throw new IllegalArgumentException("Chromosome not found in index");
        }

        ArrayList<Long> data = chromosomeData.get(chr);
        long chrStart = data.get(0);
        int charsPerLine = data.get(1).intValue();
        int charsPerLineWithNewline = data.get(2).intValue();

        // Calculate the file pointer position
        int lineNumber = (start - 1) / charsPerLine;
        int lineOffset = (start - 1) % charsPerLine;
        long filePointer = chrStart + lineNumber * charsPerLineWithNewline + lineOffset;

        // Seek to the calculated position
        raf.seek(filePointer);

        // Read the sequence
        StringBuilder sequence = new StringBuilder();
        int length = end - start + 1;
        int remaining = length;

        while (remaining > 0) {
            int toRead = Math.min(charsPerLine - lineOffset, remaining);
            byte[] buffer = new byte[toRead];
            raf.read(buffer, 0, toRead);
            sequence.append(new String(buffer));

            remaining -= toRead;
            lineOffset = 0;

            if (remaining > 0) {
                raf.skipBytes(charsPerLineWithNewline - charsPerLine);
            }
        }

        String seq = sequence.toString();
        if (!strand) {
            //seq = getReverseComplement(seq);
        }

        return seq;
    }
}