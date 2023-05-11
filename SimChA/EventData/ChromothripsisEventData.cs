// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ChromothripsisEventData : ContigEventData
{
    public List<long> StopsList { get; }
    public List<int> SelectionList { get; }
    
    public ChromothripsisEventData(Random rnd, CNEventPars cnEventPars, int contigId, long contigLen) : base(cnEventPars, contigId)
    {
        ContigId = contigId;
        double chromothripsisLen = cnEventPars.Get("Size", 100_000_000L);
        int shardCount = Sampling.GetFragCount(rnd, contigLen / chromothripsisLen);
        StopsList = Sampling.GetStopsForShards(rnd, contigLen, shardCount);
        int shardsKept = rnd.Next(1, shardCount);
        SelectionList = Enumerable.Range(0, shardCount).Shuffle(rnd).Take(shardsKept).ToList();
    }
    
    public override void ApplyEvent(Karyotype kar)
        => kar.ApplyChromothripsis(ContigId, StopsList, SelectionList);

    public override string ToString()
        => $"contig:{ContigId};" +
           $"stops:{string.Join(",", StopsList)};" +
           $"selection:{string.Join(",", SelectionList)}";
}