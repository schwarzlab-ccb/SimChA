namespace SimChA.EventData;

public record CNEventDesc(
    CNEventType EventType,
    int Depth,
    string Description = "",
    double DeltaFitness = 0,
    double TotalFitness = 0,
    double Time = 0)
{
    public static string Header() 
        => "event_type" +
           "\tdepth" +
           "\tdescription" +
           "\tdelta_fitness" +
           "\ttotal_fitness" +
           "\ttime";
    
    public string ToTSV() =>
        $"{EventType}" +
        $"\t{Depth}" +
        $"\t{Description}" +
        $"\t{DeltaFitness:f4}" +
        $"\t{TotalFitness:f4}" +
        $"\t{Time:f4}";
}