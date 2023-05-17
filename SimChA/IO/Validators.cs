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
            case CNEventType.Translocation:
            case CNEventType.ChromDeletion:
            case CNEventType.ChromDuplication:
            case CNEventType.BreakageFusionBridge:
            case CNEventType.WholeGenomeDoubling:
            case CNEventType.TailDeletion:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalDuplication:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.Chromothripsis:
            case CNEventType.Chromoplexy:
            case CNEventType.TIChain:
            case CNEventType.TIBridge:
            case CNEventType.TICycle:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
                break;
            default:
                throw new ArgumentOutOfRangeException($"Unknown event type {cnEventPars.Type}");
        }
    }
    
    // Removes signatures with probability <=0 and validates the rest
    public static List<Signature> ValidateSignatures(List<Signature>? signatures)
    {
        if (signatures is null || signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        foreach (var sig in signatures.Where(sig => sig.Prob > 0 && sig.Events.Any(e => e.Prob > 0)))
        {
            if (sig.Events is null || sig.Events.Count == 0)
            {
                throw new Exception($"Signature {sig.Id} does not have any events.");
            }
            double probSum = sig.Events.Sum(e => e.Prob);
            if (probSum <= 0)
            {
                throw new Exception($"Signature {sig.Id} has a total probability of {probSum}.");
            }
            foreach(var cnEventP in sig.Events)
            {
                ValidateEvent(cnEventP);
            }
        }
        return signatures;
    }
}