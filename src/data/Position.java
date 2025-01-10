package data;

import java.io.Serializable;

public record Position(
        double horizontal,
        double vertical
)implements Serializable {
    public Position(){
        this(0, 0);
    }
}