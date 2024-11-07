package ReadSimulator;

import java.util.HashMap;

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
}