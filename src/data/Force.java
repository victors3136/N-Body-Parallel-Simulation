package data;

public record Force (
    double horizontal,
    double vertical
){
    public Force(){
        this(0,0);
    }
}

