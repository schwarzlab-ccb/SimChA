using System.Runtime.InteropServices;
using SimChA.Computation;
using SimChA.Simulation;

namespace SimChA.EventData;

public record ArmEventData : ContigEventData
{
    public int CentromereIndex { get; }
    public bool PArm { get; }
    
    // Constructor used for Arm-level events CNEventPars
    public ArmEventData(Random rnd, CNEventPars CNEventPars, int contigId, int centromereIndex, bool pArm) : base(CNEventPars, contigId)
    {
        CentromereIndex = centromereIndex;
        PArm = pArm;
    }

    public override void ApplyEvent(Karyotype kar)
    {
        if (EventType == CNEventType.ArmDeletion)
        {
            kar.ApplyArmDeletion(ContigId, CentromereIndex, PArm);
        }
        else if (EventType == CNEventType.ArmDuplication)
        {
            kar.ApplyArmDuplication(ContigId, CentromereIndex, PArm);
        }
        else
        {
            throw new InvalidOperationException("Invalid event type for ArmEventData");
        }
    }
    
    public override string ToString()
        => PArm ? $"contig:{ContigId};p;cent:{CentromereIndex}" : $"contig:{ContigId};q;cent:{CentromereIndex}";
}
