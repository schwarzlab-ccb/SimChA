using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointMutationData : ContigEventData
{
    public long Location { get; }
    public Nucleotide NewNucleotide {get;}
    public Nucleotide OldNucleotide {get;}
    public PointMutationData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Location = Sampling.GetInternalPos(rnd, contigLen);
        // TODO: look this position up from the assembly, then later, how do I see it from the
        // previous alterations
        OldNucleotide = Nucleotide.A;
        NewNucleotide = Sampling.SampleNucleotide(rnd);
    }

    public override void ApplyEvent(Karyotype kar)
        => kar.ApplySNV(ContigId, Location, NewNucleotide);

    public override string ToString()
        => EventType switch
        {
            CNEventType.SNV => $"contig:{ContigId};location:{Location};old:{OldNucleotide},new:{NewNucleotide}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}