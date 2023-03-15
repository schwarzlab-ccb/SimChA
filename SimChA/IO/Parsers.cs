using System.Globalization;
using System.Text.Json;
using System.Text.RegularExpressions;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.IO;

public static class Parsers
{
    public static SimParams ParseSimParams(string serializedJSON)
    {
        SimParams? res;
        var options = new JsonSerializerOptions { IncludeFields = true };
        try
        {
            res = JsonSerializer.Deserialize<SimParams>(serializedJSON, options);
            if (res is null)
            {
                throw new Exception($"Could not parse the simulation parameters:\n{serializedJSON}");
            }
        }
        catch (JsonException)
        {
            throw new Exception($"Could not parse the simulation parameters:\n{serializedJSON}");
        }
        if (res.Seed < 0)
        {
            return res with { Seed = new Random().Next() };
        }
        return res;
    }

    public static void ValidateSignatures(List<Signature>? signatures)
    {
        if (signatures is null || signatures.Count == 0)
        {
            throw new Exception("No signatures were provided.");
        }
        foreach (var sig in signatures)
        {
            if (sig.Events is null || sig.Events.Count == 0)
            {
                throw new Exception($"Signature {sig.Id} does not have any events.");
            }
            foreach(var ev in sig.Events)
            {
                switch (ev.Type)
                {
                    case CNEventType.Translocation:
                    case CNEventType.ChromDeletion:
                    case CNEventType.ChromDuplication:
                    case CNEventType.BreakageFusionBridge:
                    case CNEventType.WholeGenomeDoubling:
                        break;
                    case CNEventType.TailDeletion:
                    case CNEventType.InternalDeletion:
                    case CNEventType.InternalDuplication:
                    case CNEventType.InternalInversion:
                    case CNEventType.InvertedDuplication:
                        if (ev.Params == null || !ev.Params.ContainsKey("Mean"))
                        {
                            throw new Exception($"The mean of {ev.Type} in signature {sig.Id} not specified.");
                        }
                        break;
                    case CNEventType.Chromothripsis:
                        break;
                    case CNEventType.Chromoplexy:
                        break;
                    default:
                        throw new ArgumentOutOfRangeException($"Unknown event type {ev.Type}");
                }
            }
        }
    }

    private static Karyotype MakeKaryotype(List<Region> regionsA, List<Region> regionsB, List<GenRange> missingRanges, bool sexXX)
    {
        regionsA = RegionOps.GlueNeighbours(regionsA);
        var contigA = new Contig(regionsA);
        regionsB = RegionOps.GlueNeighbours(regionsB);
        var contigB = new Contig(RegionOps.GlueNeighbours(regionsB));
        // Add full missing chromosomes
        var chrPresent = HGRef.ChrIDsForSex(sexXX).ToDictionary(c => c, c => false);
        regionsA.ForEach(r => chrPresent[r.ChrNo] = true);
        regionsB.ForEach(r => chrPresent[r.ChrNo] = true);
        var missingChrs = chrPresent.Where(pair => !pair.Value).Select(pair => new GenRange(0, HGRef.GetChromLen(pair.Key), pair.Key));
        var totalMissing = missingChrs.Concat(missingRanges).ToList();
        return new Karyotype(new List<Contig> {contigA, contigB}, totalMissing, sexXX);
    }
    
    public static Dictionary<string, Karyotype> ParseCNAProfile(TextReader cnaFile)
    {
        Dictionary<string, Karyotype> result = new();
        var missingRanges = new List<GenRange>();
        var regionsA = new List<Region>();
        var regionsB = new List<Region>();
        bool sexXX = true;
        var lastSample = "";
        string sample = "";
        var lastChr = ChrNo.chr1;
        var lastPos = 0L;
        cnaFile.ReadLine(); // Skip header
        while (cnaFile.ReadLine() is { } line)
        {
            string[] lineSplit = line.Split('\t');
            sample = lineSplit[0];

            // Set the new sample
            if (sample != lastSample)
            {
                if (regionsA.Any() || regionsB.Any())
                {
                    result[lastSample] = MakeKaryotype(regionsA, regionsB, missingRanges, sexXX);
                }
                // Reset
                regionsA.Clear();
                regionsB.Clear();
                missingRanges.Clear();
                lastSample = sample;
                lastChr = ChrNo.chr1;
                lastPos = 0L;
                sexXX = true;
            }

            // Parse the line
            var chrNo = (ChrNo) Enum.Parse(typeof(ChrNo), lineSplit[1]);
            if (chrNo == ChrNo.chrY)
            {
                sexXX = false;
            }
            int start = int.Parse(lineSplit[2]) - 1;
            int end = int.Parse(lineSplit[3]);
            int cnA = int.Parse(lineSplit[4]);
            int cnB = int.Parse(lineSplit[5]);
            
            // Check for missing ranges
            if (chrNo == lastChr)
            {
                // Range skipped
                if (lastPos != start)
                {
                    missingRanges.Add(new GenRange(lastPos, start, lastChr));
                }
            }
            else
            {
                // Till the end of a chromosome
                if (lastPos != HGRef.GetChromLen(lastChr))
                {
                    missingRanges.Add(new GenRange(lastPos, HGRef.GetChromLen(lastChr), lastChr));
                }
                // Start of a chromosome
                if (start != 0)
                {
                    missingRanges.Add(new GenRange(0, start, chrNo));
                }
            }
            lastChr = chrNo;
            lastPos = end;
            
            // Add the new regions
            for (int i = 0; i < cnA; i++)
            {
                regionsA.Add(new Region(start, end, new ChrID(chrNo, true)));
            }
            for (int i = 0; i < cnB; i++)
            {
                regionsB.Add(new Region(start, end, new ChrID(chrNo, false)));
            }
        }
        
        // Add the last sample
        result[lastSample] = MakeKaryotype(regionsA, regionsB, missingRanges, sexXX);
        return result;
    }

    public static Dictionary<ChrNo, List<Gene>> ParseGeneList(TextReader geneFile, bool isFemale)
    {
        // Pre-initialization
        var noEnum = Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToList();
        var geneList = noEnum.ToDictionary(c => c, c => new List<Gene>());
        while (geneFile.ReadLine() is { } line)
        {
            if (line == "") continue;
            string[] genString = line.Split('\t');
            //Don't include Y chromosome in genes list if clone is female
            if (isFemale && (ChrNo) Enum.Parse(typeof(ChrNo), genString[2]) == ChrNo.chrY) continue;
            string name = genString[3];
            double fitness = double.Parse(genString[4], CultureInfo.InvariantCulture.NumberFormat);
            var chrNum = (ChrNo) Enum.Parse(typeof(ChrNo), genString[0]);
            // Convert to zero-based [start, end) index 
            var region = new GenRange(int.Parse(genString[1]) - 1, int.Parse(genString[2]), chrNum);
            var gene = new Gene(name, region, fitness);
            geneList[chrNum].Add(gene);
        }
        foreach (var pair in geneList)
        {
            pair.Value.Sort((g1, g2) => g1.Range.Start.CompareTo(g2.Range.Start));
        }
        return geneList;
    }

    public static List<Clone> ParseNewick(string newickString, bool isFemale)
    {
        List<Clone> clones = new();
        if (newickString == "")
        {
            return clones;
        }
        
        const string regexPattern = @"(?<closeChildren>[(])|" +
                                    @"(?<openChildren>[)])|" +
                                    @"(?<branchLength>:+[0-9a-zA-Z-_.]+)|" +
                                    @"(?<nodeName>[0-9a-zA-Z-_.]+)|" +
                                    @"(?<nextNode>[,])|" +
                                    @"(?<root>[;])";
        //Reverse order of newick file to start with root
        const RegexOptions regexOptions = RegexOptions.RightToLeft | RegexOptions.IgnorePatternWhitespace;

        var matches = Regex.Matches(newickString, regexPattern, regexOptions);
        if (!matches.Any() || matches[0].Value != ";")
        {
            throw new Exception("Newick file is not in the right format");
        }

        bool branchLength = CheckBranchLength(matches);
        //Iterate through Regex-Matches
        var parentIds = new List<int> {-1};
        foreach (Match match in matches)
        {
            switch (match.Value)
            {
                case ";":
                    //create root node
                    var root = CreateClone(clones.Count,  parentIds.Last(), match.NextMatch(),
                        match.NextMatch().NextMatch(), isFemale, branchLength, 0);
                    clones.Add(root);
                    parentIds.Add(clones[^1].CloneId);
                    break;
                case "(":
                    //remove parent from parentList and add childrenID to parent
                    parentIds.RemoveAt(parentIds.Count - 1);
                    break;
                case ")":
                    // Add parent to parentIds and then create a child
                    parentIds.Add(clones[^1].CloneId);
                    goto case ",";
                case ",":
                    //create new child
                    int parentId = parentIds.Last();
                    var child = CreateClone(clones.Count, parentId, match.NextMatch(),
                        match.NextMatch().NextMatch(), isFemale, branchLength, clones[parentId].TotalMutations);
                    clones.Add(child);
                    clones[parentId].ChildrenIDs.Add(clones.Count - 1);
                    break;
            }
        }
        if (!clones.Any())
        {
            throw new Exception("No clones found in newick file. Might not be the right format.");
        }

        return clones;
    }

    //create clone from newick Match
    private static Clone CreateClone(int id, int parentId, Match branchLengthMatch, Match nameMatch, bool isFemale,
        bool branchLength,
        int parentMutations)
    {
        string nameClone = nameMatch.Groups["nodeName"].Value != ""
            ? nameMatch.Value
            : branchLengthMatch.Groups["nodeName"].Value != ""
                ? branchLengthMatch.Value
                : "C" + id;
        int mutCount = branchLengthMatch.Groups["branchLength"].Value != ""
            ?
            (int) Math.Ceiling(float.Parse(branchLengthMatch.Value.Remove(0, 1)))
            : branchLength
                ? 0
                : 1;
        int totalMut = parentMutations + mutCount;
        var clone = new Clone(id, parentId, nameClone, mutCount, new Karyotype(isFemale), totalMut);
        return clone;
    }

    //check for branch-length in newick file
    private static bool CheckBranchLength(MatchCollection matches)
    {
        var branchLength = false;
        foreach (Match match in matches)
        {
            if (match.Groups["branchLength"].Value != "")
            {
                branchLength = true;
                break;
            }
        }

        Console.Write(branchLength ? "" : "No branch-lengths were found, using 1 as branch-length.");
        return branchLength;
    }
}