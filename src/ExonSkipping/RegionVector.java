package ExonSkipping;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class RegionVector {
    private List<CDS> cdsList = new ArrayList<>();

    public void addCDS(CDS cds) {
        cdsList.add(cds);
        Collections.sort(cdsList, Comparator.comparingInt(CDS::getStart));
    }

    public List<CDS> getCdsList() {
        return cdsList;
    }
}