namespace SimChA.DataTypes;

public record GenomeReference(Region[] XYGenome, long XYGenomeLen, Region[] XXGenome, long XXGenomeLen);