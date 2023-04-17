// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;

namespace SimChA.EventData;

public record ChromoplexyEventData : BaseEventData
{
    public readonly List<int>? ContigIdList = new();
    public readonly List<List<long>>? StopsList = new();
    public readonly List<int>? SequenceList = new();
    public readonly List<long>? BreakpointsList = new();
    public ChromoplexyEventData(CNEventP eventP, List<int> contigIds, List<List<long>> stops, List<int> sequence, List<long> breakpoints) : base(eventP)
    {
        ContigIdList = contigIds;
        StopsList = stops;
        SequenceList = sequence;
        BreakpointsList = breakpoints;
    }

    public IEnumerable<int> GetSequence()
    {
        return SequenceList;
    }
}