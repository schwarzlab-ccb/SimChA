// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record TemplatedEventData : BaseEventData
{
    public readonly int ContigId = -1;
    public readonly List<Region>? RegionsList = new();
    public readonly int LastContigId = -1;
    public TemplatedEventData(CNEventP eventP, int contigId, List<Region> regions, int lastContigId = -1) : base(eventP)
    {
        ContigId = contigId;
        RegionsList = regions;
        LastContigId = lastContigId;
    }
    public override string ToString() => $"{EventType}\t{ContigId}\t{string.Join(",", RegionsList)}\t{LastContigId}";
}