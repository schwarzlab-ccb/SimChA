using SimChA.Computation;
using SimChA.Data;
using SimChA.Simulation;

namespace SimChA.EventData;

public record TailEventData : ContigEventData
{
    public long Length { get; }
    public bool Direction { get; }
    
    // Constructor used for Tail CNEventPars
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, long contigLen) : base(CNEventPars, contigId)
    {
        //Length = Sampling.GetPos(rnd, contigLen, CN);
        Length = Sampling.GetExpSeg(rnd, contigLen, CNEventPars.Size);
        Direction = rnd.CoinFlip();
    }
    
    public TailEventData(Random rnd, CNEventPars CNEventPars, int contigId, 
        IEnumerable<(long start, long end)> centromeres) : base(CNEventPars, contigId)
    {
        var cent =  centromeres.Shuffle(rnd).First();
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
    
    public override string ToString()
        => $"contig:{ContigId};length:{Length};dir:{Direction}";
}
