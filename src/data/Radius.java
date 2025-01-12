package data;

import java.io.Serializable;

public record Radius(double value) implements Serializable {
    public Radius() {
        this(0);
    }
}
