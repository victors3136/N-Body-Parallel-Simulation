package data;

import java.util.ArrayList;
import java.util.List;

public record Quadrant(
        int index,
        int innerQuadrantCount,
        Mass mass,
        Position bottomLeftCorner,
        Position centerOfMass,
        Dimension dimensions,
        List<Quadrant> innerQuadrants) {
    Quadrant(int index, Dimension dimensions, Position bottomLeftCorner) {
        this(index,
                0,
                new Mass(),
                bottomLeftCorner,
                new Position(),
                dimensions,
                new ArrayList<>());
    }
}
