namespace SimChA.EventData;

public record CNEventDesc(
    BaseEventData EventData, 
    int Depth, 
    double DeltaFitness = 0, 
    double TotalFitness = 0, 
    int NumRejections = 0, 
    // Instrumentation for the slot-outcome question: of the tries a slot spent, how many were lost
    // to the drawn event type being impossible on this karyotype (GenerateCNEventData returning
    // null), and how many to the proposal being inviable (an essential gene at zero copies). Both
    // currently `continue` in EvoSimulator, so their rates are invisible -- and the impossibility
    // rate is what a slot-consuming exit would turn into skips, which is the number needed to size
    // that change before making it. NumRejections stays the count of tries lost to fitness.
    int NumImpossible = 0,
    int NumInviable = 0,
    // Why the slot ended: Accepted, Impossible, Rejected (fitness exhausted the try budget) or
    // PloidyCeiling. Today a Skip record cannot say which of the last three it was.
    string Outcome = "",
    // The event type drawn on the final try. Recorded for skips too, where it is currently lost --
    // per-type skip rates are the diagnostic for whether the configured event mix survives.
    string AttemptedType = "",
    string Signature = "", 
    string RegionsGained = "",
    string RegionsLost = "", 
    string Karyotype = "",
    int? ContigSlotsBefore = null,
    int? ContigSlotsAfter = null)
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
           "\tnum_impossible" +
           "\tnum_inviable" +
           "\toutcome" +
           "\tattempted_type" +
           "\tsignature" + 
            (PrintDelta ? "\tregions_gained\tregions_lost" : "") +
            (PrintKaryotype ? "\tkaryotype" : "") +
            "\tcontig_slots_before\tcontig_slots_after";

    public string ToTSV() =>
        $"{EventData.EventType}" +
        $"\t{Depth}" +
        $"\t{EventData.EventDesc()}" +
        $"\t{DeltaFitness:f4}" +
        $"\t{TotalFitness:f4}" +
        $"\t{NumRejections}" +
        $"\t{NumImpossible}" +
        $"\t{NumInviable}" +
        $"\t{Outcome}" +
        $"\t{AttemptedType}" +
        $"\t{Signature}" +
        (PrintDelta ? $"\t{RegionsGained}\t{RegionsLost}" : "") +
        (PrintKaryotype ? $"\t{Karyotype}" : "") +
        $"\t{ContigSlotsBefore}\t{ContigSlotsAfter}";
}