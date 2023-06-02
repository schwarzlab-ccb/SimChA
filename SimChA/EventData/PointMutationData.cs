using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointMutationData : ContigEventData
{
    public long Location { get; }
    public Nucleotide MutatedNucleotide {get;}
    public PointMutationData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Location = Sampling.GetInternalPos(rnd, contigLen);
        MutatedNucleotide = Sampling.SampleNucleotide(rnd);
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.SNV)
        {
            kar.ApplySNV(ContigId, Location, MutatedNucleotide);
        }
        else if (EventType == CNEventType.PointInsertion)
        {
            kar.ApplyPointInsertion(ContigId, Location, MutatedNucleotide);
        }
        else if (EventType == CNEventType.PointDeletion)
        {
            kar.ApplyPointDeletion(ContigId, Location);
        }
    }
    public override string ToString()
        => EventType switch
        {
            CNEventType.PointInsertion or CNEventType.SNV => $"contig:{ContigId};location:{Location};inserted:{MutatedNucleotide}",
            CNEventType.PointDeletion  => $"contig:{ContigId};location:{Location}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}