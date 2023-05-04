// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record PairEventData : BaseEventData
{
    public readonly int ContigIdA = -1;
    public readonly int ContigIdB = -1;
    public readonly long PosA = -1;
    public readonly long PosB = -1;
    public readonly bool Direction;
    
    // Constructor used for Translocation
    public PairEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigA, int contigB) : base(eventP)
    {
        ContigIdA = contigA;
        long lenA = kar.ContigLen(contigA);
        PosA = Sampling.GetInternalPos(rnd, lenA);
        
        ContigIdB = contigB;
        long lenB = kar.ContigLen(contigB);
        PosB = Sampling.GetInternalPos(rnd, lenB);
        
        double invProb = eventP.Get("InvProb", 0.0);
        Direction = invProb != 0.0 && rnd.CoinFlip(invProb);
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyTranslocation(ContigIdA, ContigIdB, PosA, PosB, Direction);
    
    public override string ToString()
        => $"contigA:{ContigIdA};contigB:{ContigIdB};posA:{PosA};posB:{PosB};dir:{Direction}";
}