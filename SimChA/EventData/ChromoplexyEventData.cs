// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ChromoplexyEventData : BaseEventData
{
    public List<int> ContigIds { get; } // Contigs that are involved in the event
    public List<List<long>> Stops { get; } // Breakpoints for each contig
    public List<int> Sequence { get; } // Sequence of shards (given by the number of pieces in each contig)
    public List<long> Breakpoints { get; } // Breakpoints for the whole event - shards are first joined then split again
    
    // TODO: Validate
    public ChromoplexyEventData(Random rnd, CNEventP eventP, List<(int id, long len)> seq) : base(eventP)
    {
        int contigCount = Math.Min(seq.Count, Sampling.GetChromoplexySiteCount(rnd));
        double size = eventP.Get("Size", 10_000_000L);
        ContigIds = new List<int>();
        Stops = new List<List<long>>(); 
        var totalLen = 0L;
        var totalFrags = 0;
        for (int i = 0; i <contigCount; i++)
        {
            ContigIds.Add(seq[i].id);
            totalLen += seq[i].len;
            int partsCount = Sampling.GetFragCount(rnd, seq[i].len / size);
            totalFrags += partsCount;
            Stops.Add(Sampling.GetStopsForShards(rnd, seq[i].len, partsCount));
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