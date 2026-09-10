namespace SimChA.IO;

public record EvoParams(
    int MaxTries = 100, // Maximum number of tries to generate a new set of events in a time step
    double Decay = 0.25 // Decay factor for fitness matching mode: higher values make acceptance stricter as fewer events remain
);
