namespace SimChA.Data;

public record SNV(long Location, string ChrNo, Nucleotide Alt)
{
    public override string ToString() => $"chr:{ChrNo};location:{Location};alt:{Alt}";
}