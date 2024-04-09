// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.DataTypes;
using SimChA.EventData;

namespace SimChA.IO;

public static class Validators
{
    private static void ValidateEvent(CNEventPars cnEventPars)
    {
        switch (cnEventPars.Type)
        {
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
            case CNEventType.WholeGenomeDoubling:
            case CNEventType.BreakageFusionBridge:
            case CNEventType.TailDeletion:
            case CNEventType.TailDuplication:
            case CNEventType.SNV:
            case CNEventType.ArmDeletion:
            case CNEventType.ArmDuplication:
                break;
                
            case CNEventType.InternalDeletion:
            case CNEventType.InternalDuplication:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.Translocation:
            case CNEventType.Chromothripsis:
            case CNEventType.Chromoplexy:
                if (cnEventPars.Size <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Size\": 1000000");
                break;
            
            case CNEventType.TIChain:
            case CNEventType.TIBridge:
            case CNEventType.TICycle:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
                if (cnEventPars.Size <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Size\": 1000000");
                if (cnEventPars.Frag <= 0) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Frag\": 5");
                break;
            
            default:
                throw new ArgumentOutOfRangeException($"Unknown event type {cnEventPars.Type}");
        }
    }
    
    // Removes signatures with probability <=0 and validates the rest
    public static void ValidateSignatures(Dictionary<string, Signature> signatures)
    {
        foreach (var sig in signatures.Where(sig => sig.Value.Prob > 0 && sig.Value.Events.Any(e => e.Prob > 0)))
        {
            if (sig.Value.Events is null || sig.Value.Events.Count == 0)
            {
                throw new Exception($"Signature {sig.Key} does not have any events.");
            }
            double probSum = sig.Value.Events.Sum(e => e.Prob);
            if (probSum <= 0)
            {
                throw new Exception($"Signature {sig.Key} has a total probability of {probSum}.");
            }
            foreach(var cnEventP in sig.Value.Events)
            {
                ValidateEvent(cnEventP);
            }
        }
    }
}
