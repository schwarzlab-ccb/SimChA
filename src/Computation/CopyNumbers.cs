using SimChA.Data;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static Dictionary<string, List<int>> GetJointSegmentation(List<string> chromNames, List<Sample> samples)
    {
        var breaks = samples.Select(s => s.Karyotype.CalcBreaks()).ToList();
        var segmentation = chromNames.ToDictionary(
            chrom => chrom, 
            chrom => breaks.SelectMany(br => br.GetValueOrDefault(chrom, [])).ToHashSet().OrderBy(b => b).ToList());
        return segmentation;
    }

    public static List<CopyNumber> CalcCNs(IDictionary<string, List<int>> allBreaks, List<Contig> contigs)
    {
        var result = new List<CopyNumber>();
        foreach ((string chrom, var breaks) in allBreaks)
        {
            for (int i = 0; i < breaks.Count - 1; i++)
            {
                int start = breaks[i];
                int end = breaks[i + 1];
                var seg = new GenRange(start, end, chrom);
                var cns = contigs.Select(c => c.GetCNs(seg));
                (int cnA, int cnB, int nSNVs) = cns.Aggregate((CNA: 0, CNB: 0, SNV: 0), (acc, vals)
                    => (acc.CNA + vals.CNA, acc.CNB + vals.CNB, acc.SNV + vals.SNV));
                var cn = new CopyNumber(start, end, chrom, cnA, cnB, nSNVs);
                result.Add(cn);
            }
        }
        return result;
    }

    public static List<CopyNumber> CalcCNs(Karyotype karyotype, IDictionary<string, List<int>>? breaks = null)
        => karyotype.CalcCNs(breaks ?? karyotype.CalcBreaks());

    private static void AddEndpoint(
        Dictionary<string, Dictionary<long, (int H1, int H2)>> endpoints,
        string chrom,
        long position,
        int deltaH1,
        int deltaH2)
    {
        if (!endpoints.TryGetValue(chrom, out var chromEndpoints))
        {
            chromEndpoints = [];
            endpoints.Add(chrom, chromEndpoints);
        }

        var current = chromEndpoints.GetValueOrDefault(position);
        var updated = (H1: current.H1 + deltaH1, H2: current.H2 + deltaH2);
        if (updated == (0, 0))
        {
            chromEndpoints.Remove(position);
        }
        else
        {
            chromEndpoints[position] = updated;
        }
    }

    private static void AddDeltaInterval(
        List<CopyNumber> deltas,
        ref CopyNumber? run,
        string chrom,
        long start,
        long end,
        int deltaH1,
        int deltaH2)
    {
        if (deltaH1 == 0 && deltaH2 == 0)
        {
            return;
        }
        if (run is not null && run.End == start && run.CNH1 == deltaH1 && run.CNH2 == deltaH2)
        {
            run = new CopyNumber(run.Start, end, chrom, deltaH1, deltaH2, 0);
            return;
        }
        if (run is not null)
        {
            deltas.Add(run);
        }
        run = new CopyNumber(start, end, chrom, deltaH1, deltaH2, 0);
    }

    // Each region represents one H1 or H2 copy over [start, end), so its coverage can be
    // encoded as +coverage at start and -coverage at end. Karyotype and Contig traverse
    // their private children and pass only those scalar endpoints here. Contributions from
    // before use -1 and contributions from after use +1; therefore, the running sum between
    // sorted endpoints is exactly the after-minus-before copy-number delta. Equal adjacent
    // intervals are merged and zero intervals are omitted. This avoids rebuilding two complete
    // copy-number tables (and their unused SNV counts) after every simulated event.
    public static List<CopyNumber> DiffKaryotypes(Karyotype before, Karyotype after)
    {
        var endpoints =
            new Dictionary<string, Dictionary<long, (int H1, int H2)>>();
        void AddEndpointToDelta(string chrom, long position, int deltaH1, int deltaH2)
            => AddEndpoint(endpoints, chrom, position, deltaH1, deltaH2);
        before.AccumulateCopyNumberEndpoints(AddEndpointToDelta, -1);
        after.AccumulateCopyNumberEndpoints(AddEndpointToDelta, 1);

        var deltas = new List<CopyNumber>();
        foreach (string chrom in before.ChromNames
                     .Concat(after.ChromNames)
                     .Distinct())
        {
            if (!endpoints.TryGetValue(chrom, out var chromEndpoints) || chromEndpoints.Count == 0)
            {
                continue;
            }

            long? previous = null;
            int deltaH1 = 0;
            int deltaH2 = 0;
            CopyNumber? run = null;
            foreach (var endpoint in chromEndpoints.OrderBy(pair => pair.Key))
            {
                if (previous is long start && start < endpoint.Key)
                {
                    AddDeltaInterval(
                        deltas,
                        ref run,
                        chrom,
                        start,
                        endpoint.Key,
                        deltaH1,
                        deltaH2);
                }
                deltaH1 += endpoint.Value.H1;
                deltaH2 += endpoint.Value.H2;
                previous = endpoint.Key;
            }
            if (run is not null)
            {
                deltas.Add(run);
            }
        }
        return deltas;
    }

    public static double CalcPloidy(RefGen refGen, List<CopyNumber> copyNumbers, SexType sex)
        => 2 * copyNumbers.Select(c => c.Length * (c.CNH1 + c.CNH2)).Sum() / (float) refGen.SexGenomeLen[(int)sex];
}
