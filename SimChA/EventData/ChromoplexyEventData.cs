// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using Extreme.Mathematics;
using SimChA.DataTypes;
using SimChA.Misc;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ChromoplexyEventData : BaseEventData
{
    public List<int> ContigIdList { get; }
    public List<List<long>> StopsList { get; }
    public List<int> SequenceList { get; }
    public List<long> BreakpointsList{ get; }
    
    public ChromoplexyEventData(Random rnd, Karyotype kar, CNEventP eventP, IEnumerator<int> IDsEnumerator) : base(eventP)
    {
        int chrCount = Sampling.GetChromoplexySiteCount(rnd);
        ContigIdList = new List<int>();
        StopsList = new List<List<long>>();
        long contigLen = kar.ContigLen(IDsEnumerator.Current); 
        var totalLen = 0L;
        var totalFrags = 0;
        for (var i = 0; i < chrCount; i++, IDsEnumerator.MoveNext())
        {
            ContigIdList.Add(IDsEnumerator.Current);
            long thisLen = kar.ContigLen(IDsEnumerator.Current); 
            totalLen += thisLen;
            int partsCount = Sampling.GetFragCount(rnd, thisLen / (double) contigLen);
            totalFrags += partsCount;
            StopsList.Add(Sampling.GetStopsForShards(rnd, contigLen, partsCount));
        }
        SequenceList = Enumerable.Range(0, totalFrags).Shuffle(rnd).ToList();
        BreakpointsList = Sampling.GetStopsForShards(rnd, totalLen, chrCount);
    }

    public override string ApplyEvent(Karyotype kar)
        => kar.ApplyEvent(this);
}