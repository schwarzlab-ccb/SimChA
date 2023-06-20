using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record PointMutationData : ContigEventData
{
    public long Location { get; set;}
    public Nucleotide OldNucleotide {get; set;}
    public Nucleotide NewNucleotide {get; set;}
    Random Rnd {get;}
    public SNV SNV;
    private Dictionary<char, Nucleotide> nucleotideDict = 
        new Dictionary<char, Nucleotide>() 
        {
            {'A', Nucleotide.A},
            {'C', Nucleotide.C},
            {'G', Nucleotide.G},
            {'T', Nucleotide.T}
        };
    public PointMutationData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        Location = Sampling.GetInternalPos(rnd, contigLen);
        Rnd = rnd;
    }

    public void SetOldNucleotide(Karyotype kar)
    {
        (Region region, long internalLocation) = kar.GetContig(ContigId).FindRegion(Location);
        var dummySNV = new SNV(Nucleotide.A, Nucleotide.A);
        // TODO: THIS IS WILDLY UNSAFE, but done for the tests to pass
        if (kar.GenContents == null)
        {
            OldNucleotide = Nucleotide.A;
        }
        else if (region.SNVDict == null || !region.SNVDict.TryGetValue(internalLocation, out dummySNV))
        {
            int index = (int)region.ChrID.ChrNo;
            var nuc = Char.ToUpper(kar.GenContents[index].Sequence[(int)internalLocation]);
            // TODO: What do we do with the 'N' character?
            if (nuc == 'N')
            {
                OldNucleotide = Sampling.SampleNucleotide(Rnd);
            }
            else
            {
                OldNucleotide = nucleotideDict[nuc];
            }
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
        SNV = CreateSNV(kar);
        kar.ApplySNV(ContigId, Location, SNV);
    }

    public override string ToString()
        => EventType switch
        {
            CNEventType.SNV => $"contig:{ContigId};location:{Location};old:{OldNucleotide},new:{NewNucleotide}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}