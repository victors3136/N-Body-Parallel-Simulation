package data;

public class Point {
    private Position position;
    private final double radius;
    private final Mass mass;
    private Velocity velocity;
    private Force force;

    public Point() {
        this(new Position(), 0, new Mass(), new Velocity(), new Force());
    }

    public Point(Position position, double radius, Mass mass, Velocity velocity) {
        this(position, radius, mass, velocity, new Force());
    }

    public Point(Position position, double radius, Mass mass, Velocity velocity, Force force) {
        this.position = position;
        this.radius = radius;
        this.mass = mass;
        this.velocity = velocity;
        this.force = force;
    }

    public Position position() {
        return position;
    }

    public double radius() {
        return radius;
    }

    public Mass mass() {
        return mass;
    }

    public Velocity velocity() {
        return velocity;
    }

    public Force force() {
        return force;
    }

    public void setForce(Force force) {
        this.force = force;
    }

    public void setVelocity(Velocity velocity) {
        this.velocity = velocity;
    }

    public void setPosition(Position position) {
        this.position = position;
    }
}
