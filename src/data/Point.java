package data;

import java.io.Serializable;

public class Point implements Serializable {
    private Position position;
    private Radius radius;
    private Mass mass;
    private Velocity velocity;
    private Force force;

    public Point() {
        this(new Position(), new Radius(), new Mass(), new Velocity(), new Force());
    }

    public Point(Position position, Radius radius, Mass mass, Velocity velocity) {
        this(position, radius, mass, velocity, new Force());
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

    public void setMass(Mass mass) {
        this.mass = mass;
    }
    public void setRadius(Radius radius) {
        this.radius = radius;
    }

    public Force directForce(Point other) {
        final var dx = other.position.horizontal() - this.position.horizontal();
        final var dy = other.position.vertical() - this.position.vertical();
        final var distance = this.position.distance(other.position);

        final var magnitude = (
                CommonCore.GravitationalConstant * mass.value() * other.mass.value())
                / (distance * distance);
        final var factor = magnitude / distance;
        final var fx = factor * dx;
        final var fy = factor * dy;
        return new Force(fx, fy);
    }

    public void updatePosition() {
        final var fx = 0;
        final var fy = 0;

        final var vx = velocity.horizontal() + force.horizontal() / mass.value() * CommonCore.timeIncrement;
        final var vy = velocity().vertical() + force.vertical() / mass.value() * CommonCore.timeIncrement;

        final var px = position.horizontal() + velocity.horizontal() * CommonCore.timeIncrement;
        final var py = position.vertical() + velocity.vertical() * CommonCore.timeIncrement;

        force = new Force(fx, fy);
        velocity = new Velocity(vx, vy);
        position = new Position(px, py);
    }


}

