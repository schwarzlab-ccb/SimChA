namespace SimChA.DataTypes;

public record CopyNumber(Region Segment, int CNH1, int CNH2, int NSNVs)
{
    public string ToTSV()
        => string.Join('\t', Segment.ChrID.ChrNo, Segment.Start, Segment.End, CNH1, CNH2, NSNVs);
}