package data;

import java.io.Serializable;

public record Dimension (
        double horizontal,
        double vertical
) implements Serializable {
    Dimension(){
        this(0, 0);
    }
}
