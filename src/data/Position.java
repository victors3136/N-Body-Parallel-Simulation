package data;

import java.io.Serializable;

public record Position(
        double horizontal,
        double vertical
) implements Serializable {
    public Position(){
        this(0, 0);
    }

    public double distance(Position that) {
        return Math.sqrt(
                Math.pow(this.horizontal - that.horizontal, 2)
                        + Math.pow(this.vertical - that.vertical, 2));
    }
}