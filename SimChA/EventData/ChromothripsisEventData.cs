// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record ChromothripsisEventData : BaseEventData
{
    public readonly int ContigId = -1;
    public readonly List<long>? StopsList = new();
    public readonly List<int>? SelectionList = new();
    public ChromothripsisEventData(CNEventP eventP, int contigId, List<long> stops, List<int> selection) : base(eventP)
    {
        ContigId = contigId;
        StopsList = stops;
        SelectionList = selection;
    }

    public IEnumerable<int> GetSelection()
    {
        return SelectionList;
    }
}