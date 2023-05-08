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
    public readonly bool Inverted;
    
    // Constructor used for Translocation
    public PairEventData(Random rnd, CNEventP eventP, int contigA, long lenA, int contigB, long lenB) : base(eventP)
    {
        ContigIdA = contigA;
        PosA = Sampling.GetInternalPos(rnd, lenA);
        ContigIdB = contigB;
        PosB = Sampling.GetInternalPos(rnd, lenB);
        double invProb = eventP.Get("InvProb", 0.0);
        Inverted = invProb != 0.0 && rnd.CoinFlip(invProb);
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyTranslocation(ContigIdA, ContigIdB, PosA, PosB, Inverted);
    
    public override string ToString()
        => $"contigA:{ContigIdA};contigB:{ContigIdB};posA:{PosA};posB:{PosB};dir:{Inverted}";
}