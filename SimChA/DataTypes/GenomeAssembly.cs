using System.Text.Json.Serialization;

namespace SimChA.DataTypes;

[JsonConverter(typeof(JsonStringEnumConverter))]
public enum GenomeAssembly
{
    hg19,
    hg38
}