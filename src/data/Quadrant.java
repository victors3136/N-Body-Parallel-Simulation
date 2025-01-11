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
    private Point innerPoint;

    public Quadrant() {
        this(0, 0, new Mass(), new Position(), new Position(), new Dimension(), new ArrayList<>());
    }

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
        this.innerPoint = null;
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

    public void insert(Point point) {
        if (innerQuadrantCount == 0 && innerPoint == null) {
            innerPoint = point;
        } else if (innerQuadrantCount == 0) {
            subdivide();
            assert innerPoint != null;
            final var toMove = innerPoint;
            subquadrantContaining(toMove).insert(toMove);
            subquadrantContaining(point).insert(point);

        } else {
            subquadrantContaining(point).insert(point);
        }
        updateMass(point);
    }

    private void updateMass(Point point) {
        final var newMass = mass.value() + point.mass().value();
        final var cX = (centerOfMass.horizontal() * mass.value() + point.position().horizontal() * point.mass().value()) / newMass;
        final var cY = (centerOfMass.vertical() * mass.value() + point.position().vertical() * point.mass().value()) / newMass;
        final var newCenterOfMass = new Position(cX, cY);
        mass = new Mass(newMass);
        centerOfMass = newCenterOfMass;
    }

    private Quadrant subquadrantContaining(Point point) {
        final var anchorX = bottomLeftCorner.horizontal();
        final var anchorY = bottomLeftCorner.vertical();
        final var halfHeight = dimensions.vertical() / 2;
        final var halfWidth = dimensions.horizontal() / 2;
        final var position = point.position();
        final var x = position.horizontal();
        final var y = position.vertical();
        if (x < anchorX || y < anchorY) {
            throw new RuntimeException("How did we get here?");
        }
        if ((x < anchorX + halfWidth) && (y < anchorY + halfHeight)) {
            return innerQuadrants.get(RelativePosition.BOTTOM_LEFT.ordinal());
        }
        if ((x <= anchorX + halfWidth) && (y <= anchorY + 2 * halfHeight)) {
            return innerQuadrants.get(RelativePosition.TOP_LEFT.ordinal());
        }
        if ((x <= anchorX + 2 * halfWidth) && (y <= anchorY + 2 * halfHeight)) {
            return innerQuadrants.get(RelativePosition.TOP_RIGHT.ordinal());
        }
        if ((x <= anchorX + 2 * halfWidth) && (y <= anchorY + halfHeight)) {
            return innerQuadrants.get(RelativePosition.BOTTOM_RIGHT.ordinal());
        }
        throw new RuntimeException("How did we get here?");
    }

    private void subdivide() {
        final var anchorX = bottomLeftCorner.horizontal();
        final var anchorY = bottomLeftCorner.vertical();
        final var halfHeight = dimensions.vertical() / 2;
        final var halfWidth = dimensions.horizontal() / 2;
        {
            Quadrant bottomLeft = new Quadrant();
            bottomLeft.setBottomLeftCorner(new Position(anchorX, anchorY));
            bottomLeft.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants.add(bottomLeft);
        }
        {
            Quadrant topLeft = new Quadrant();
            topLeft.setBottomLeftCorner(new Position(anchorX, anchorY + halfHeight));
            topLeft.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants.add(topLeft);
        }
        {
            Quadrant topRight = new Quadrant();
            topRight.setBottomLeftCorner(new Position(anchorX + halfWidth, anchorY + halfHeight));
            topRight.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants.add(topRight);
        }
        {
            Quadrant bottomRight = new Quadrant();
            bottomRight.setBottomLeftCorner(new Position(anchorX + halfWidth, anchorY));
            bottomRight.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants.add(bottomRight);
        }
        innerQuadrantCount = 4;
    }
}
