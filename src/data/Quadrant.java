package data;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class Quadrant implements Serializable {
    private int index;
    private int innerQuadrantCount;
    private Mass mass;
    private Position bottomLeftCorner;
    private Position centerOfMass;
    private Dimension dimensions;
    private List<Quadrant> innerQuadrants;

    public Quadrant(int index, Dimension dimensions, Position bottomLeftCorner) {
        this(index,
                0,
                new Mass(),
                bottomLeftCorner,
                new Position(),
                dimensions,
                new ArrayList<>());
    }

    public Quadrant(int index,
             int innerQuadrantCount,
             Mass mass,
             Position bottomLeftCorner,
             Position centerOfMass,
             Dimension dimensions,
             List<Quadrant> innerQuadrants) {
        this.index = index;
        this.innerQuadrantCount = innerQuadrantCount;
        this.mass = mass;
        this.bottomLeftCorner = bottomLeftCorner;
        this.centerOfMass = centerOfMass;
        this.dimensions = dimensions;
        this.innerQuadrants = innerQuadrants;
    }

    public int index() {
        return index;
    }

    public int innerQuadrantCount() {
        return innerQuadrantCount;
    }

    public Mass mass() {
        return mass;
    }

    public Position bottomLeftCorner() {
        return bottomLeftCorner;
    }

    public Position centerOfMass() {
        return centerOfMass;
    }

    public Dimension dimensions() {
        return dimensions;
    }

    public List<Quadrant> innerQuadrants() {
        return innerQuadrants;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public void setInnerQuadrantCount(int innerQuadrantCount) {
        this.innerQuadrantCount = innerQuadrantCount;
    }

    public void setMass(Mass mass) {
        this.mass = mass;
    }

    public void setBottomLeftCorner(Position bottomLeftCorner) {
        this.bottomLeftCorner = bottomLeftCorner;
    }

    public void setCenterOfMass(Position centerOfMass) {
        this.centerOfMass = centerOfMass;
    }

    public void setDimensions(Dimension dimensions) {
        this.dimensions = dimensions;
    }

    public void setInnerQuadrants(List<Quadrant> innerQuadrants) {
        this.innerQuadrants = innerQuadrants;
    }
}
