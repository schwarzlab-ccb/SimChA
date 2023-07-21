namespace SimChA.DataTypes;

public record CopyNumber(GenRange Segment, int CNH1, int CNH2, int NSNVs)
{
    public string ToTSV()
        => string.Join('\t', Segment.ChrNo, Segment.Start + 1, Segment.End, CNH1 >= 0 ? CNH1 : "NA", CNH2 >= 0 ? CNH2 : "NA", NSNVs >= 0 ? NSNVs : "NA");
}
