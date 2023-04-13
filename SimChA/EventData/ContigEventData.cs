// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.EventData;

public record ContigEventData : BaseEventData
{
    public readonly int ContigId = -1;
    
    // Constructor used for whole-chromosome events
    public ContigEventData(CNEventP eventP, int contigId)  : base(eventP)
    {
        ContigId = contigId;
    }
    
    public override string ToString() => $"{EventType}\t{ContigId}";
    
}