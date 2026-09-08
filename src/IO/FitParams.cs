// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Computation;

namespace SimChA.IO;

public record FitParams
(
    double Stress = 0,
    double TsgOg = 0,
    double Essentiality = 0,
    // Reject any event that would leave an essential gene at zero copies, instead of pricing the
    // loss. Losing both copies is inviability rather than unfitness, and it lets Essentiality mean
    // one thing -- the cost of haploinsufficiency -- rather than setting the balance between two
    // states. Only EvoSimulator enforces it; MonteCarlo has no acceptance step to reject at and
    // FitnessMatching optimises a different target. Set false to price the loss instead, which is
    // what every fit before this did. Configs predating the field deserialise to this default.
    bool ProhibitEssentialLoss = true,
    string GeneSet = "Empty"
);
