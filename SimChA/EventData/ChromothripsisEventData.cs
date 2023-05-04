// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ChromothripsisEventData : BaseEventData
{
    public int ContigId { get; }
    public List<long> StopsList { get; }
    public List<int> SelectionList { get; }
    
    public ChromothripsisEventData(Random rnd, Karyotype kar, CNEventP eventP, int contigId) : base(eventP)
    {
        ContigId = contigId;
        long contigLen = kar.ContigLen(ContigId);
        double chromothripsisLen = eventP.Get("Size", 100_000_000L);
        int shardCount = Sampling.GetFragCount(rnd, contigLen / chromothripsisLen);
        StopsList = Sampling.GetStopsForShards(rnd, contigLen, shardCount);
        int shardsKept = rnd.Next(1, StopsList.Count);
        SelectionList = Enumerable.Range(0, shardCount).Shuffle(rnd).Take(shardsKept).ToList();
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyChromothripsis(ContigId, StopsList, SelectionList);

    public override string ToString()
        => $"contig:{ContigId};" +
           $"stops:{string.Join(",", StopsList)};" +
           $"selection:{string.Join(",", SelectionList)}";
}