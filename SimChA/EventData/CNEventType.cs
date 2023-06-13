// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Text.Json.Serialization;

namespace SimChA.EventData;

[JsonConverter(typeof(JsonStringEnumConverter))]
public enum CNEventType
{
    ChromDeletion,
    ChromDuplication,
    TailDeletion,
    InternalDeletion,
    InternalDuplication,
    InternalInversion,
    InvertedDuplication,
    Translocation,
    BreakageFusionBridge,
    WholeGenomeDoubling,
    Chromothripsis,
    Chromoplexy,
    TIChain,
    TICycle,
    TIBridge,
    Pyrgo,
    Rigma,
    SNV
    // Tyfonas
}