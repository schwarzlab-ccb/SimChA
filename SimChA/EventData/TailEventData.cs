using SimChA.Computation;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public bool Direction { get; }
    
    // Constructor used for Tail CNEventPars
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) 
        : base(CNEventPars, contigId, contigLen)
    {
        //Length = Sampling.GetPos(rnd, contigLen, CN);
        Length = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Frac);
        Direction = rnd.CoinFlip();
    }
    
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, IEnumerable<(long start, long end)> cents)
        : base(CNEventPars, contigId, -1)
    {
        var cent =  cents.Shuffle(rnd).First();
        Length = rnd.CoinFlip() ? cent.start : cent.end;
        Direction = rnd.CoinFlip();
    }

    public override void ApplyEvent(Karyotype kar)
    {
        switch (EventType)
        {
            case CNEventType.TailDeletion:
            case CNEventType.ArmDeletion:
                kar.ApplyTailDeletion(ContigId, Length, Direction);
                break;
            case CNEventType.TailDuplication:
            case CNEventType.ArmDuplication:
                kar.ApplyTailDuplication(ContigId, Length, Direction);
                break;
            case CNEventType.BreakageFusionBridge:
                kar.ApplyBFB(ContigId, Length, Direction);
                break;
            default:
                throw new Exception($"Invalid event type {EventType} for TailEventData");
        }
    }
    
    public override string EventDesc()
        => base.EventDesc() + "loc:" + (Direction ? "start;" : "end;");
}
