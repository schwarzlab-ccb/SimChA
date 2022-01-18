using SimChA.Computation;

namespace SimChA.DataTypes;
using MathNet.Numerics.Distributions;

public static class RawData
{
    public static List<SNPData> CalcSingleSubclone(
        List<CopyNumber> copyNumbers, List<SNP> snps, float purity = 1f, float baferror = 0.1F,
        float readstd = 2f, int readdepth = 10)
    {
        // TODO: implement purity
        var result = new List<SNPData>();
        int curSegmentId = 0;
        foreach (var snp in snps)
        {
            curSegmentId = copyNumbers.FindIndex(curSegmentId, cn =>
                ReferenceGenome.ChromosomeStartMap[cn.Segment.ChromId.ChromNum] + cn.Segment.End > snp.AbsPos);

            if (copyNumbers[curSegmentId].CNH1 + copyNumbers[curSegmentId].CNH2 <= 0)
                continue;

            int readsH1 = (int)Math.Round(Math.Max(0, Normal.Sample(copyNumbers[curSegmentId].CNH1 * readdepth, readstd)));
            int readsH2 = (int)Math.Round(Math.Max(0, Normal.Sample(copyNumbers[curSegmentId].CNH2 * readdepth, readstd)));

            float logr = (float)Math.Log2((readsH1 + readsH2) / (2f * readdepth));
            float baf = snp.Heterozygous && readsH1 + readsH2 > 0 ? (float)readsH2 / (readsH1 + readsH2) : -1;

            result.Add(new SNPData(snp, logr, baf));

            // TODO: Extra noise for logr and baf
        }
        return result;
    }

    public static void CalcMultipleSubclones()
    {
        // TODO: Should take multiple subclones as input and mix them
        throw new NotImplementedException();
    }

    private static string FirstLine(bool isFirst, string type)
        => isFirst ? $"\tchrom\tpos\t{type}\n" : "";

    public static string PrintBAF(List<SNPData> rawData)
        => FirstLine(true, "BAF") + string.Join("\n", rawData.Select(r => r.PrintBAF()));

    public static string PrintLogR(List<SNPData> rawData)
        => FirstLine(true, "logR") + string.Join("\n", rawData.Select(r => r.PrintLogR()));
}