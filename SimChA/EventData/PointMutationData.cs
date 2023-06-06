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
        if (EventType != CNEventType.PointDeletion)
        {
            if (EventType == CNEventType.SNV)
            {
                OldNucleotide = Nucleotide.A;
                // TODO: how do I look this up from the assembly, and then later,
                // how do I see it from the previous alterations
            }
            NewNucleotide = Sampling.SampleNucleotide(rnd);
        }        
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.SNV)
        {
            kar.ApplySNV(ContigId, Location, NewNucleotide);
        }
        else if (EventType == CNEventType.PointInsertion)
        {
            kar.ApplyPointInsertion(ContigId, Location, NewNucleotide);
        }
        else if (EventType == CNEventType.PointDeletion)
        {
            kar.ApplyPointDeletion(ContigId, Location);
        }
    }
    public override string ToString()
        => EventType switch
        {
            CNEventType.SNV => $"contig:{ContigId};location:{Location};old:{OldNucleotide},new:{NewNucleotide}",
            CNEventType.PointInsertion => $"contig:{ContigId};location:{Location};inserted:{NewNucleotide}",
            CNEventType.PointDeletion  => $"contig:{ContigId};location:{Location}",
            _ => throw new ArgumentOutOfRangeException(nameof(EventType), EventType, null)
        };
        
}