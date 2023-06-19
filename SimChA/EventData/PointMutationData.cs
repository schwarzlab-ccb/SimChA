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
        var region = kar.GetContig(ContigId).FindRegion(Location);
        var dummySNV = new SNV(Nucleotide.A, Nucleotide.C);
        if (region.SNVDict == null || !region.SNVDict.TryGetValue(Location, out dummySNV))
        {
            OldNucleotide = Nucleotide.A;
        }
        else
        {
            OldNucleotide = dummySNV.NewNucleotide;
        }
    }

    public SNV CreateSNV(Karyotype kar)
    {
        SetOldNucleotide(kar);
        NewNucleotide = Sampling.SampleNucleotide(Rnd, OldNucleotide);
        return new SNV(OldNucleotide, NewNucleotide);
    }
    public override void ApplyEvent(Karyotype kar)
    {
        var SNV = CreateSNV(kar);
        kar.ApplySNV(ContigId, Location, SNV);
    }

    public override string ToString()
        => EventType switch
        {
            CNEventType.SNV => $"contig:{ContigId};location:{Location};old:{OldNucleotide},new:{NewNucleotide}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}