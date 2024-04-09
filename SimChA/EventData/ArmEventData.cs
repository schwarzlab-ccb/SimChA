using System.Runtime.InteropServices;
using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ArmEventData : ContigEventData
{
    public bool PArm { get; }
    
    // Constructor used for Tail CNEventPars
    public ArmEventData(Random rnd, CNEventPars CNEventPars, int contigId) : base(CNEventPars, contigId)
    {
        PArm = rnd.CoinFlip();
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.ArmDeletion)
        {
            kar.ApplyArmDeletion(ContigId, PArm);
        }
        else if (EventType == CNEventType.ArmDuplication)
        {
            kar.ApplyArmDuplication(ContigId, PArm);
        }
        else
        {
            throw new InvalidOperationException("Invalid event type for ArmEventData");
        }
    }
    
    public override string ToString()
        => PArm ? $"contig:{ContigId};p" : $"contig:{ContigId};q";
}
