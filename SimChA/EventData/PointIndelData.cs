using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointIndelData : ContigEventData
{
    public long Location { get; }
    public Nucleotide InsertedNucleotide {get;}

    // Constructor used for point insertions and deletions
    public PointIndelData(Random rnd, CNEventPars CNEventPars, int contigId, long location) : base(CNEventPars, contigId)
    {
        Location = location;
        InsertedNucleotide = 
            (EventType == CNEventType.PointInsertion) ? Sampling.SampleNucleotide(rnd) : Nucleotide.A ;
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.PointInsertion)
        {
            kar.ApplyPointInsertion(ContigId, Location, InsertedNucleotide);
        }
        else if (EventType == CNEventType.PointDeletion)
        {
            kar.ApplyPointDeletion(ContigId, Location);
        }
    }

    public override string ToString()
        => EventType switch
        {
            CNEventType.PointInsertion => $"contig:{ContigId};location:{Location};inserted:{InsertedNucleotide}",
            CNEventType.PointDeletion  => $"contig:{ContigId};location:{Location}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}