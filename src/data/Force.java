package data;

public record Force (
    double horizontal,
    double vertical
){
    Force(){
        this(0,0);
    }
}

