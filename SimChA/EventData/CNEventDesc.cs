namespace SimChA.EventData;

public record CNEventDesc(BaseEventData EventData, int Depth, double DeltaFitness = 0, double TotalFitness = 0, int NumTries = 0)
{
    public static string Header()
        => "event_type" +
           "\tdepth" +
           "\tdescription" +
           "\tdelta_fitness" +
           "\ttotal_fitness" +
           "\tnum_tries";

    public string ToTSV() =>
        $"{EventData.EventType}" +
        $"\t{Depth}" +
        $"\t{EventData.EventDesc()}" +
        $"\t{DeltaFitness:f4}" +
        $"\t{TotalFitness:f4}" +
        $"\t{NumTries}";
}