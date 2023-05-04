// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ChromoplexyEventData : BaseEventData
{
    public List<int> ContigIds { get; }
    public List<List<long>> Stops { get; }
    public List<int> Sequence { get; }
    public List<long> Breakpoints { get; }
    
    // TODO: Validate
    public ChromoplexyEventData(Random rnd, Karyotype kar, CNEventP eventP, IEnumerable<int> ids) : base(eventP)
    {
        int contigCount = Sampling.GetChromoplexySiteCount(rnd);
        double size = eventP.Get("Size", 10_000_000L);
        ContigIds = ids.Take(contigCount).ToList();
        Stops = new List<List<long>>(); 
        var totalLen = 0L;
        var totalFrags = 0;
        foreach (int id in ContigIds)
        {
            long thisLen = kar.ContigLen(id); 
            totalLen += thisLen;
            int partsCount = Sampling.GetFragCount(rnd, thisLen / size);
            totalFrags += partsCount;
            Stops.Add(Sampling.GetStopsForShards(rnd, id, partsCount));
        }
        Sequence = Enumerable.Range(0, totalFrags).Shuffle(rnd).ToList();
        Breakpoints = Sampling.GetStopsForShards(rnd, totalLen, contigCount);
    }

    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyChromoplexy(ContigIds, Stops, Sequence, Breakpoints);

    public override string ToString()
        => $"contigs:[{string.Join(",", ContigIds)}];" +
           $"stops:[{string.Join(",", Stops.Select(s => $"[{string.Join(",", s)}]"))}];" +
           $"sequence:[{string.Join(",", Sequence)}];" +
           $"breakpoints:[{string.Join(",", Breakpoints)}]";
}