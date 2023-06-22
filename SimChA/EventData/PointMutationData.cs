using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointMutationData : ContigEventData
{
    public long Location { get; }
    public Nucleotide Base { get; }

    public PointMutationData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) 
        : base(CNEventPars, contigId)
    {
        Location = Sampling.GetInternalPos(rnd, contigLen);
        Base = Sampling.SampleBase(rnd);
    }
    
    public override void ApplyEvent(Karyotype kar)
    {
        kar.ApplySNV(ContigId, Location, Base);
    }

    public override string ToString()
        => $"contig:{ContigId};location:{Location};base:{Base}";
}