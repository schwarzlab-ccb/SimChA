namespace SimChA.EventData;

public record CNEventDesc(
    BaseEventData EventData, 
    int Depth, 
    double DeltaFitness = 0, 
    double TotalFitness = 0, 
    int NumRejections = 0, 
    string Signature = "", 
    string RegionsGained = "",
    string RegionsLost = "", 
    string Karyotype = "")
{
    public static bool PrintDelta { get; set; }
    public static bool PrintKaryotype  { get; set; }
    
    public static string Header()
        => "event_type" +
           "\tdepth" +
           "\tdescription" +
           "\tdelta_fitness" +
           "\ttotal_fitness" +
           "\tnum_rejections" +
           "\tsignature" + 
            (PrintDelta ? "\tregions_gained\tregions_lost" : "") +
            (PrintKaryotype ? "\tkaryotype" : "");

    public string ToTSV() =>
        $"{EventData.EventType}" +
        $"\t{Depth}" +
        $"\t{EventData.EventDesc()}" +
        $"\t{DeltaFitness:f4}" +
        $"\t{TotalFitness:f4}" +
        $"\t{NumRejections}" +
        $"\t{Signature}" +
        (PrintDelta ? $"\t{RegionsGained}\t{RegionsLost}" : "") +
        (PrintKaryotype ? $"\t{Karyotype}" : "");
}