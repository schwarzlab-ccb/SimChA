// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Simulation;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record SNVData : ContigEventData
{
    public long Location { get; }
    public Nucleotide MutatedNucleotide {get;}
    public SNVData(Random rnd, CNEventPars CNEventPars, int contigId, long location, Nucleotide mutatedNucleotide) : base(CNEventPars, contigId)
    {
        Location = location;
        MutatedNucleotide = mutatedNucleotide;
    }
    public override void ApplyEvent(Karyotype kar)
    {
        kar.ApplySNV(ContigId, Location, MutatedNucleotide);
    }
    
    public override string ToString() 
        => $"contig:{ContigId};location:{Location};";
}