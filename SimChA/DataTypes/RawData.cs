using SimChA.Computation;

namespace SimChA.DataTypes;
using MathNet.Numerics.Distributions;

public static class RawData
{
    public static List<SNPData> CalcSingleSubclone(
        List<CopyNumber> copyNumbers, List<SNP> snps, bool isfemale, float purity = 1f, float baferror = 0.1F,
        float readstd = 2f, int readdepth = 10, float gamma = 1f)
    {
        var result = new List<SNPData>();
        int curSegmentId = 0;
        float tumorPloidy = CopyNumbers.CalcPloidy(copyNumbers, isfemale);
        foreach (var snp in snps)
        {
            curSegmentId = copyNumbers.FindIndex(curSegmentId, cn =>
                ReferenceGenome.ChromosomeStartMap[cn.Segment.ChromId.ChromNum] + cn.Segment.End > snp.AbsPos);

            if (copyNumbers[curSegmentId].CNH1 + copyNumbers[curSegmentId].CNH2 <= 0)
                continue;

            int readsNormal = (int)Math.Round(Math.Max(0, Normal.Sample(2 * readdepth, readstd)));
            int readsH1 = (int)Math.Round(Math.Max(0,
                Normal.Sample(copyNumbers[curSegmentId].CNH1 * readdepth, readstd)));
            int readsH2 = (int)Math.Round(Math.Max(0,
                Normal.Sample(copyNumbers[curSegmentId].CNH2 * readdepth, readstd)));

            float logr = gamma * (float)Math.Log2((1 - purity) * readdepth * 2 + purity * (readsH1 + readsH2) / (readdepth * 2 * (1 - purity) + readdepth * purity * tumorPloidy));
            float baf = snp.Heterozygous && readsH1 + readsH2 > 0 ? (float)(((1 - purity) * readsNormal + purity * readsH2) / ((1 - purity) * 2 * readsNormal + purity * (readsH1 + readsH2))) : -1;

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