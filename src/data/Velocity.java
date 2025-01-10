package data;

public record Velocity(
        double horizontal,
        double vertical
) {
    public Velocity() {
        this(0, 0);
    }
}
