using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointMutationData : ContigEventData
{
    public long Location { get; }
    public Nucleotide OldNucleotide {get; set;}
    public Nucleotide NewNucleotide {get; set;}
    Random Rnd {get;}
    public PointMutationData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Location = Sampling.GetInternalPos(rnd, contigLen);
        Rnd = rnd;
    }

    public void SetOldNucleotide(Karyotype kar)
    {
        OldNucleotide = Nucleotide.A;
    }
    public SNV CreateSNV()
    {
        // TODO: look this position up from the assembly, then later, how do I see it from the
        // previous alterations
        NewNucleotide = Sampling.SampleNucleotide(Rnd, OldNucleotide);
        return new SNV(OldNucleotide, NewNucleotide);
    }

    // ApplyEvent for SNVs has to be a little more involved because there's no way of setting the 
    public override void ApplyEvent(Karyotype kar)
    {
        SetOldNucleotide(kar);
        var SNV = CreateSNV();
        kar.ApplySNV(ContigId, Location, SNV);
    }

    public override string ToString()
        => EventType switch
        {
            CNEventType.SNV => $"contig:{ContigId};location:{Location};old:{OldNucleotide},new:{NewNucleotide}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}