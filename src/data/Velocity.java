package data;

import java.io.Serializable;

public record Velocity (
        double horizontal,
        double vertical
)implements Serializable {
    public Velocity() {
        this(0, 0);
    }
}
