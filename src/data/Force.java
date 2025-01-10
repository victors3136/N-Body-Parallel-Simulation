package data;

import java.io.Serializable;

public record Force (
    double horizontal,
    double vertical
)implements Serializable {
    public Force(){
        this(0,0);
    }
}

