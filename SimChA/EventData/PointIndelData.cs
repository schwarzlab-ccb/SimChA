using SimChA.Simulation;

namespace SimChA.EventData;

public record PointIndelData : ContigEventData
{
    public long Location { get; }

    // Constructor used for point insertions and deletions
    public PointIndelData(Random rnd, CNEventPars CNEventPars, int contigId, long location) : base(CNEventPars, contigId)
    {
        Location = location;
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.PointInsertion)
        {
            kar.ApplyPointInsertion(ContigId, Location);
        }
        else if (EventType == CNEventType.PointDeletion)
        {
            kar.ApplyPointDeletion(ContigId, Location);
        }
    }

    public override string ToString()
        => $"contig:{ContigId};location:{Location}";
}