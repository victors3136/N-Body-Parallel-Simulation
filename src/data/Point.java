package data;

public class Point {
    private final Position position;
    private final double radius;
    private final Mass mass;
    private final Velocity velocity;
    private Force force;

    public Point(Position position, double radius, Mass mass, Velocity velocity) {
        this.position = position;
        this.radius = radius;
        this.mass = mass;
        this.velocity = velocity;
        force = new Force();
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
}
