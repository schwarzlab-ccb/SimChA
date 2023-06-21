using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointMutationData : ContigEventData
{
    public long Location { get; set;}
    
    public SNV SNV { get; private set; }
    
    public PointMutationData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Location = Sampling.GetInternalPos(rnd, contigLen);
        var newNucleotide = Sampling.SampleNucleotide(rnd);
        SNV = new SNV(Nucleotide.N, newNucleotide);
    }

    private Nucleotide GetOldNucleotide(Karyotype kar)
    {
        var oldNuc = Nucleotide.N;
        if (kar.GenContents != null)
        {
            (var region, long internalLocation) = kar.GetContig(ContigId).FindRegion(Location);
            int index = (int) region.ChrID.ChrNo;
            char nuc = char.ToUpper(kar.GenContents[index].Sequence[(int) internalLocation]);
            oldNuc = Enum.Parse<Nucleotide>(nuc.ToString());
        }
        return oldNuc;
    }
    
    public override void ApplyEvent(Karyotype kar)
    {
        var oldNuc = GetOldNucleotide(kar);
        SNV = SNV with {OldNucleotide = oldNuc};
        kar.ApplySNV(ContigId, Location, SNV);
    }

    public override string ToString()
        => EventType switch
        {
            CNEventType.SNV => $"contig:{ContigId};location:{Location};old:{SNV.OldNucleotide},new:{SNV.NewNucleotide}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}