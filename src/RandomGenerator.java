import data.*;

import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

class RandomGenerator {
    static private final Random random;

    static {
        random = new Random();
    }

    public static Position nextPosition() {

        final var innerHorizontalBound = CommonCore.maxWidth - CommonCore.maxBodyRadius;
        final var innerVerticalBound = CommonCore.maxHeight - CommonCore.maxBodyRadius;
        return new Position(
                random.nextDouble() % innerHorizontalBound,
                random.nextDouble() % innerVerticalBound);
    }

    public static Force nextForce() {
        return new Force(
                random.nextDouble() % (CommonCore.maxInitialForce) + CommonCore.minInitialForce,
                random.nextDouble() % (CommonCore.maxInitialForce) + CommonCore.minInitialForce);
    }

    public static Mass nextMass() {
        return new Mass(random.nextDouble() % CommonCore.maxMass);
    }

    public static Radius nextRadius() {
        return new Radius(random.nextDouble() % CommonCore.maxBodyRadius);
    }

    public static Velocity nextVelocity() {
        return new Velocity(
                random.nextDouble() % CommonCore.maxInitialVelocity,
                random.nextDouble() % CommonCore.maxInitialVelocity
        );
    }

    public static Point nextPoint() {
        return new Point(nextPosition(), nextRadius(), nextMass(), nextVelocity(), nextForce());
    }


    public static List<Position> positions(int size) {
        assert size > 0;
        return IntStream.range(0, size).mapToObj(_ -> nextPosition()).toList();
    }

    public static List<Force> forces(int size) {
        assert size > 0;
        return IntStream.range(0, size).mapToObj(_ -> nextForce()).toList();
    }

    public static List<Mass> masses(int size) {
        assert size > 0;
        return IntStream.range(0, size).mapToObj(_ -> nextMass()).toList();
    }

    public static List<Velocity> velocities(int size) {
        assert size > 0;
        return IntStream.range(0, size).mapToObj(_ -> nextVelocity()).toList();
    }

    public static List<Point> points(int size) {
        assert size > 0;
        return IntStream.range(0, size).mapToObj(_ -> nextPoint()).toList();
    }
}
