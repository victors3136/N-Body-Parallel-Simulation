package data;

import java.io.Serializable;

public class Point implements Serializable {
    private Position position;
    private Radius radius;
    private Mass mass;
    private Velocity velocity;
    private Force force;

    public Point(Position position, Radius radius, Mass mass, Velocity velocity) {
        this.position = position;
        this.radius = radius;
        this.mass = mass;
        this.velocity = velocity;
        force = new Force();
    }

    public Point(Position position, Radius radius, Mass mass, Velocity velocity, Force force) {
        this.position = position;
        this.radius = radius;
        this.mass = mass;
        this.velocity = velocity;
        this.force = force;
    }

    public Position position() {
        return position;
    }

    public Radius radius() {
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
