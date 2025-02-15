package data;

import java.io.Serializable;

public class Quadrant implements Serializable {
    private int innerQuadrantCount;
    private Mass mass;
    private Position bottomLeftCorner;
    private Position centerOfMass;
    private Dimension dimensions;
    private final Quadrant[] innerQuadrants;
    private Point innerPoint;

    private static final int BOTTOM_LEFT = 0;
    private static final int TOP_LEFT = 1;
    private static final int TOP_RIGHT = 2;
    private static final int BOTTOM_RIGHT = 3;

    public Quadrant() {
        this(0, new Mass(), new Position(), new Position(), new Dimension(), new Quadrant[4]);
    }

    public Quadrant(Dimension dimensions, Position bottomLeftCorner) {
        this(
                0,
                new Mass(),
                bottomLeftCorner,
                new Position(),
                dimensions,
                new Quadrant[4]);
    }

    public Quadrant(int innerQuadrantCount,
                    Mass mass,
                    Position bottomLeftCorner,
                    Position centerOfMass,
                    Dimension dimensions,
                    Quadrant[] innerQuadrants) {
        this.innerQuadrantCount = innerQuadrantCount;
        this.mass = mass;
        this.bottomLeftCorner = bottomLeftCorner;
        this.centerOfMass = centerOfMass;
        this.dimensions = dimensions;
        this.innerQuadrants = innerQuadrants;
        this.innerPoint = null;
    }

    public Mass mass() {
        return mass;
    }

    public void setMass(Mass mass) {
        this.mass = mass;
    }

    public void setBottomLeftCorner(Position bottomLeftCorner) {
        this.bottomLeftCorner = bottomLeftCorner;
    }

    public void setDimensions(Dimension dimensions) {
        this.dimensions = dimensions;
    }

    public synchronized void insert(Point point) {
        Quadrant current = this;

        while (true) {
            if (current.innerQuadrantCount == 0 && current.innerPoint == null) {
                current.innerPoint = point;
                current.updateMass(point);
                return;
            } else if (current.innerQuadrantCount == 0) {
                current.subdivide();
                final var toMove = current.innerPoint;
                current.innerPoint = null;
                Quadrant targetQuadrantForExisting = current.subquadrantContaining(toMove);
                targetQuadrantForExisting.innerPoint = toMove;
                targetQuadrantForExisting.updateMass(toMove);
            }

            Quadrant targetQuadrantForNew = current.subquadrantContaining(point);
            if (targetQuadrantForNew.innerPoint == null && targetQuadrantForNew.innerQuadrantCount == 0) {
                targetQuadrantForNew.innerPoint = point;
                targetQuadrantForNew.updateMass(point);
                return;
            }

            current = targetQuadrantForNew;
        }
    }


    private synchronized void updateMass(Point point) {
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
        final var midX = anchorX + halfWidth;
        final var midY = anchorY + halfHeight;
        if ((x <= midX) && (y <= midY)) {
            return innerQuadrants[BOTTOM_LEFT];
        } else if ((x <= midX) && (y >= midY)) {
            return innerQuadrants[TOP_LEFT];
        } else if (y <= midY) {
            return innerQuadrants[BOTTOM_RIGHT];
        } else {
            return innerQuadrants[TOP_RIGHT];
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
            innerQuadrants[BOTTOM_LEFT] = bottomLeft;
        }
        {
            Quadrant topLeft = new Quadrant();
            topLeft.setBottomLeftCorner(new Position(anchorX, anchorY + halfHeight));
            topLeft.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants[TOP_LEFT] = topLeft;
        }
        {
            Quadrant topRight = new Quadrant();
            topRight.setBottomLeftCorner(new Position(anchorX + halfWidth, anchorY + halfHeight));
            topRight.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants[TOP_RIGHT] = topRight;
        }
        {
            Quadrant bottomRight = new Quadrant();
            bottomRight.setBottomLeftCorner(new Position(anchorX + halfWidth, anchorY));
            bottomRight.setDimensions(new Dimension(halfWidth, halfHeight));
            innerQuadrants[BOTTOM_RIGHT] = bottomRight;
        }
        innerQuadrantCount = 4;
    }

    public void addForceActingOn(Point point) {
        synchronized (point) {
            var fx = point.force().horizontal();
            var fy = point.force().vertical();

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
    }

    public Point asPoint() {
        Point virtualPoint = new Point();
        virtualPoint.setPosition(centerOfMass);
        virtualPoint.setMass(mass);
        return virtualPoint;
    }
}

