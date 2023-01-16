using SimChA.DataTypes;
using static Extreme.Statistics.Distributions.NormalDistribution;

namespace SimChA.Computation;

public static class SNPMetrics
{
    public static List<SNPData> CalcSingleSubClone(
        Random rnd,
        List<CopyNumber> copyNumbers,
        List<SNP> snps,
        bool isfemale,
        float purity = 1f,
        float baferror = 0.1f,
        float readstd = 2f,
        int readdepth = 10,
        float gamma = 1f)
    {
        var result = new List<SNPData>();
        int curSegmentId = 0;
        float tumorPloidy = CopyNumbers.CalcPloidy(copyNumbers, isfemale);
        foreach (var snp in snps)
        {
            curSegmentId = copyNumbers.FindIndex(curSegmentId, cn =>
                ReferenceGenome.ChromosomeStartMap[cn.Segment.ChrID.ChrNo] + cn.Segment.End > snp.AbsPos);

            if (copyNumbers[curSegmentId].CNH1 + copyNumbers[curSegmentId].CNH2 <= 0)
            {
                continue;
            }

            int readsNormal = (int)Math.Round(Math.Max(0,
                Sample(rnd, 2 * readdepth, readstd)));
            int readsH1 = (int)Math.Round(Math.Max(0,
                Sample(rnd, copyNumbers[curSegmentId].CNH1 * readdepth, readstd)));
            int readsH2 = (int)Math.Round(Math.Max(0,
                Sample(rnd, copyNumbers[curSegmentId].CNH2 * readdepth, readstd)));

            float compPurity = 1f - purity;
            int readsTotal = readsH1 + readsH2;
            float r = (2 * compPurity * readdepth + purity * readsTotal) /
                      (2 * compPurity * readdepth + readdepth * purity * tumorPloidy);
            float logr = gamma * (float)Math.Log2(r);
            float baf = snp.Heterozygous && readsTotal > 0
                ? (compPurity * readsNormal + purity * readsH2) /
                  (2 * compPurity * readsNormal + purity * readsTotal)
                : -1;

            result.Add(new SNPData(snp, logr, baf));
            // TODO: Extra noise for logr and baf
        }

        return result;
    }

    public static void CalcMultipleSubclones()
    {
        // TODO: Should take multiple clones as input and mix them
        throw new NotImplementedException();
    }

    private static string Header(string type)
        => $"\tchr\tpos\t{type}\n";

    public static string PrintBAF(List<SNPData> rawData)
        => Header("BAF") + string.Join("\n", rawData.Select(r => r.PrintBAF()));

    public static string PrintLogR(List<SNPData> rawData)
        => Header("logR") + string.Join("\n", rawData.Select(r => r.PrintLogR()));
}