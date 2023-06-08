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
                break;
                
            case CNEventType.BreakageFusionBridge:
            case CNEventType.TailDeletion:
            case CNEventType.InternalDeletion:
            case CNEventType.InternalDuplication:
            case CNEventType.InternalInversion:
            case CNEventType.InvertedDuplication:
            case CNEventType.Chromothripsis:
                if (cnEventPars.Pars == null || !cnEventPars.Pars.ContainsKey("Size")) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Pars\": {{\"Size\": 1000000}}");
                break;

            case CNEventType.Translocation:
                if (cnEventPars.Pars == null || !cnEventPars.Pars.ContainsKey("Size")) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Pars\": {{\"Size\": 1000000}}");
                if (cnEventPars.Pars == null || !cnEventPars.Pars.ContainsKey("Size")) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a PIvn parameter. E.g. \"Pars\": {{\"PIvn\": 0.5}}");
                break;
            
            case CNEventType.Chromoplexy:
            case CNEventType.TIChain:
            case CNEventType.TIBridge:
            case CNEventType.TICycle:
            case CNEventType.Pyrgo:
            case CNEventType.Rigma:
                if (cnEventPars.Pars == null || !cnEventPars.Pars.ContainsKey("Size")) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Size parameter. E.g. \"Pars\": {{\"Size\": 1000000}}");
                if (cnEventPars.Pars == null || !cnEventPars.Pars.ContainsKey("Frag")) 
                    throw new Exception($"Event {cnEventPars.Type} does not have a Frag parameter. E.g. \"Pars\": {{\"Frag\": 10}}");
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