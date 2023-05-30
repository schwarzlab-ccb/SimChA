namespace SimChA.DataTypes;

public record CopyNumber(GenRange Segment, int CNH1, int CNH2)
{
    public string ToTSV()
        => string.Join('\t', Segment.ChrNo, Segment.Start, Segment.End, CNH1, CNH2);
}