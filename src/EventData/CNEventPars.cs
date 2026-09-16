using System.Text.Json.Serialization;
using SimChA.Computation;

namespace SimChA.EventData;

// Shape is the second parameter of whichever length distribution the event type uses: the shape of
// the bounded Pareto for Internal* events, the alpha of the Beta for Telomere* ones. Frac fixes the
// mean of that distribution and Shape fixes its spread, so the two together determine it -- before
// Shape existed both were compile-time constants and only the mean could be configured.
//
// Null means the family default (Sampling.FixedParetoShape, Sampling.FixedBetaAlpha), which is what
// every config written before the field existed carries, so those keep their previous behaviour.
[Serializable]
public record CNEventPars(
    CNEventType Type,
    double Prob,
    double Frac = 0,
    double Frag = 0,
    string Signature = "",
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    double? Shape = null) : IHasProb;
