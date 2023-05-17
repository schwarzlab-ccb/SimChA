using System.Text.Json.Serialization;

namespace SimChA.DataTypes;

[JsonConverter(typeof(JsonStringEnumConverter))]
public enum GenomeAssembly
{
    none,
    hg19,
    hg38
}