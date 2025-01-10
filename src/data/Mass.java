package data;

import java.io.Serializable;

public record Mass (double value)implements Serializable {
    public Mass(){
        this(0);
    }
}
