// src/ExonSkipping/Intron.java
package ExonSkipping;

import java.util.Objects;

public class Intron {
    private int start;
    private int end;

    public Intron(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Intron intron = (Intron) o;
        return start == intron.start && end == intron.end;
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end);
    }
}