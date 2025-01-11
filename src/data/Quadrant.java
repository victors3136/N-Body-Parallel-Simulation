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

    public Quadrant(Dimension dimensions, Position bottomLeftCorner) {
        this(0,
                0,
                new Mass(),
                bottomLeftCorner,
                new Position(),
                dimensions,
                new ArrayList<>());
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
        final var BOTTOM_LEFT = 0;
        final var TOP_LEFT = 1;
        final var TOP_RIGHT = 2;
        final var BOTTOM_RIGHT = 3;

        final var anchorX = bottomLeftCorner.horizontal();
        final var anchorY = bottomLeftCorner.vertical();
        final var halfHeight = dimensions.vertical() / 2;
        final var halfWidth = dimensions.horizontal() / 2;
        final var position = point.position();
        final var x = position.horizontal();
        final var y = position.vertical();
        if (((x < anchorX) || (y < anchorY)) || ((x > anchorX + 2 * halfWidth) || (y > anchorY + 2 * halfHeight))) {
            throw new RuntimeException("How did we get here?");
        }
        if ((x < anchorX + halfWidth) && (y < anchorY + halfHeight)) {
            //noinspection SequencedCollectionMethodCanBeUsed
            return innerQuadrants.get(BOTTOM_LEFT);
        } else if ((x <= anchorX + halfWidth) && (y <= anchorY + 2 * halfHeight)) {
            return innerQuadrants.get(TOP_LEFT);
        } else if ((x <= anchorX + 2 * halfWidth) && (y <= anchorY + 2 * halfHeight)) {
            return innerQuadrants.get(TOP_RIGHT);
        } else {
            return innerQuadrants.get(BOTTOM_RIGHT);
        }
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

    public void addForceActingOn(Point point) {
        double fx = point.force().horizontal();
        double fy = point.force().vertical();

        if (innerQuadrantCount == 0) {
            if (innerPoint != null && !innerPoint.equals(point)) {
                final var f = point.directForce(innerPoint);
                fx += f.horizontal();
                fy += f.vertical();
                point.setForce(new Force(fx, fy));
            }
            return;
        }
        final var size = dimensions.horizontal();
        final var distance = point.position().distance(centerOfMass);
        if (size / distance < CommonCore.angle) {

            final var f = point.directForce(this.asPoint());
            fx += f.horizontal();
            fy += f.vertical();
            point.setForce(new Force(fx, fy));
            return;
        }
        for (final var innerQ : innerQuadrants) {
            if (innerQ.innerQuadrantCount != 0 || innerQ.innerPoint != null) {
                innerQ.addForceActingOn(point);
            }
        }
    }

    public Point asPoint() {
        Point virtualPoint = new Point();
        virtualPoint.setPosition(centerOfMass);
        virtualPoint.setMass(mass);
        return virtualPoint;
    }
}

