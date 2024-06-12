using System.Collections;
using System.Globalization;
using System.Text.Json;
using SimChA.DataTypes;
using SimChA.Simulation;
using System.Text;
using System.Text.RegularExpressions;

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

    // Expected format is the binned linear genome:
    // and there is a header and the columns contain:
    // SampleID, Chr, Start, End, CN hap1, CN hap2
    // Simpler version of the parser that does not calculate missing regions and does have to create karyotypes
    public static Dictionary<string, List<CopyNumber>> ParseCNAProfile(TextReader cnaFile)
    {
        Dictionary<string, List<CopyNumber>> result = new();
        cnaFile.ReadLine(); // Skip header
        while (cnaFile.ReadLine() is {} line)
        {
            string[] lineSplit = line.Split('\t');
            string sample = lineSplit[0];

            // Parse the line
            string chrNo = lineSplit[1];
            int start = int.Parse(lineSplit[2]) - 1;
            int end = int.Parse(lineSplit[3]);
            bool majorCNFound = int.TryParse(lineSplit[4], out int majorCN);
            
            bool minorCNFound = int.TryParse(lineSplit[4], out int minorCN);
            var cn = new CopyNumber(new GenRange(start, end, chrNo), 
                                    majorCNFound ? majorCN : -1, 
                                    minorCNFound ? minorCN : -1, 0);
            
            if (!result.ContainsKey(sample))
            {
                result[sample] = new List<CopyNumber>();
            }
            result[sample].Add(cn);
        }
        return result;
    }

    // Expected format is that there is a header and the columns contain:
    // SampleID, Chr, Start, End, CN hap1, CN hap2
    // NOTE: This became quite unwieldy due to the missing regions calculation,
    // however it works so don't refactor unless needed
    public static Dictionary<string, Karyotype> ParseCNAProfile(GenRef genRef, TextReader cnaFile)
    {
        Dictionary<string, Karyotype> result = new();
        var missingRanges = new List<GenRange>();
        var regionsA = new List<Region>();
        var regionsB = new List<Region>();
        var present = genRef.IncludeSexChromosomes
                    ? genRef.AllChrs.ToDictionary(c => c, _ => false)
                    : genRef.ChrIDsForAutosomes().ToDictionary(c => c, _ => false);
        string lastSample = "";
        string lastChr = genRef.AllChrs.First();
        long lastPos = 0L;
        cnaFile.ReadLine(); // Skip header
        while (cnaFile.ReadLine() is { } line)
        {
            string[] lineSplit = line.Split('\t');
            string sample = lineSplit[0];
            // Set the new sample
            if (sample != lastSample)
            {
                // First is empty
                if (regionsA.Count != 0 || regionsB.Count != 0)
                {
                    // Till the end of a chromosome
                    if (lastPos != genRef.ChrLengths[lastChr])
                    {
                        missingRanges.Add(new GenRange(lastPos, genRef.ChrLengths[lastChr], lastChr));
                    }
                    missingRanges.AddRange(
                        present
                            .Where(pair => !pair.Value)
                            .Select(c => new GenRange(0, genRef.ChrLengths[c.Key], c.Key)));
                    // Consider missing to be haplotypes by default
                    foreach (var range in missingRanges)
                    {
                        regionsA.Add(new Region(range.Start, range.End, range.ChrNo, true));
                        regionsB.Add(new Region(range.Start, range.End, range.ChrNo, false));
                    }
                    var thisKar = genRef.IncludeSexChromosomes 
                            ? new Karyotype(new List<Contig> { new(regionsA),  new(regionsB) }, missingRanges, !present[genRef.YChrName])
                            : new Karyotype(new List<Contig> { new(regionsA),  new(regionsB) }, missingRanges, false);
                    result[lastSample] = thisKar;
                }
                // Reset
                regionsA.Clear();
                regionsB.Clear();
                missingRanges.Clear();
                lastSample = sample;
                lastChr = genRef.AllChrs.First();
                lastPos = 0L;
                foreach (var pair in present)
                {
                    present[pair.Key] = false;
                }
            }
            try
            {
                // Parse the line
                string chrNo = lineSplit[1];
                present[chrNo] = true;
                int start = int.Parse(lineSplit[2]) - 1;
                int end = int.Parse(lineSplit[3]);
                int majorCN = int.Parse(lineSplit[4]);
                int minorCN = int.Parse(lineSplit[5]);

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
                    if (lastPos != genRef.ChrLengths[lastChr])
                    {
                        missingRanges.Add(new GenRange(lastPos, genRef.ChrLengths[lastChr], lastChr));
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
                for (int i = 0; i < majorCN; i++)
                {
                    regionsA.Add(new Region(start, end, chrNo, true));
                }
                for (int i = 0; i < minorCN; i++)
                {
                    regionsB.Add(new Region(start, end, chrNo, false));
                }
            }
            catch (Exception e)
            {
                throw new Exception($"Could not parse the CNA profile:\n{line}\n{e.Message}");
            }
        }

        // Till the end of a chromosome
        if (lastPos != genRef.ChrLengths[lastChr])
        {
            missingRanges.Add(new GenRange(lastPos, genRef.ChrLengths[lastChr], lastChr));
        }
        missingRanges.AddRange(
            present
                .Where(pair => !pair.Value)
                .Select(c => new GenRange(0, genRef.ChrLengths[lastChr], c.Key)));
        // Consider missing to be haplotypes by default
        foreach (var range in missingRanges)
        {
            regionsA.Add(new Region(range.Start, range.End, range.ChrNo, true));
            regionsB.Add(new Region(range.Start, range.End, range.ChrNo, false));
        }

        // Add the last sample
        var newRegs = new List<Contig> { new(regionsA), new(regionsB) };
        var kar = genRef.IncludeSexChromosomes 
                ? new Karyotype(newRegs, missingRanges, !present[genRef.YChrName])
                : new Karyotype(newRegs, missingRanges, false);
        result[lastSample] = kar;
        return result;
    }

    public static Dictionary<string, List<Gene>> ParseGeneList(TextReader geneFile, List<string> chrNames)
    {
        // Pre-initialization
        var geneList = chrNames.ToDictionary(c => c, _ => new List<Gene>());
        while (geneFile.ReadLine() is { } line)
        {
            if (line == "") continue;
            string[] genString = line.Split('\t');
            string name = genString[3];
            double fitness = double.Parse(genString[4], CultureInfo.InvariantCulture.NumberFormat);
            string chrNum = genString[0];
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

    public static List<CloneIn> ParseClones(TextReader cloneStream, bool parseFitness, string sep)
    {
        const string idKey = "ID";
        const string parentIDKey = "ParentID";
        const string distanceKey = "Distance";
        const string fitnessKey = "Fitness";

        string? firstLine = cloneStream.ReadLine() ?? throw new Exception("CloneIn file is empty.");
        var header = firstLine.Split(sep).Select(s => s.Trim()).ToList();
        var columns = new Dictionary<string, int> { { idKey, -1 }, { parentIDKey, -1 }, { distanceKey, -1 } };
        if (parseFitness)
        {
            columns.Add(fitnessKey, -1);
        }

        foreach (var column in columns)
        {
            int idx = header.IndexOf(column.Key);
            if (idx == -1) throw new Exception($"CloneIn file does not contain {column.Key} column.");
            columns[column.Key] = idx;
        }

        var cloneFitness = new List<CloneIn>();
        while (cloneStream.ReadLine() is { } line)
        {
            var lineSplit = line.Split(sep).Select(s => s.Trim()).ToList();
            int id = int.Parse(lineSplit[columns[idKey]]);
            int parentId = int.Parse(lineSplit[columns[parentIDKey]]);
            int distance = int.Parse(lineSplit[columns[distanceKey]]);
            double fitness = parseFitness
                ? double.Parse(lineSplit[columns[fitnessKey]], CultureInfo.InvariantCulture.NumberFormat)
                : 0.0;
            var clone = new CloneIn(id, parentId, distance, fitness);
            cloneFitness.Add(clone);
        }
        return cloneFitness;
    }

    public static List<(double fitness, int eventCount)> ParseClones(TextReader fitnessStream, FitnessParams fParams)
    {
        var output = new List<(double fitness, int eventCount)>();
        string? firstLine = fitnessStream.ReadLine() ?? throw new Exception("Fitness file is empty.");
        // Continue past the header
        while (fitnessStream.ReadLine() is { } line)
        {
            var lineSplit = line.Split("\t").Select(s => s.Trim()).ToList();
            double stressTerm = double.Parse(lineSplit[4], CultureInfo.InvariantCulture.NumberFormat)*fParams.Stress;
            var tsg = double.Parse(lineSplit[5], CultureInfo.InvariantCulture.NumberFormat);
            var og  = double.Parse(lineSplit[6], CultureInfo.InvariantCulture.NumberFormat);
            double tsgogTerm = (og + tsg) * fParams.TsgOg;
            double essTerm = double.Parse(lineSplit[7], CultureInfo.InvariantCulture.NumberFormat)*fParams.Essentiality;
            double totalFitness = 1.0 + (stressTerm + tsgogTerm + essTerm)*fParams.TotalStrength;
            // If the file includes data on how many chromosomal events the sample underwent
            int eventCount = lineSplit.Count >= 9 ? int.Parse(lineSplit[8]) : -1;
            output.Add((totalFitness, eventCount));
        }
        return output;
    }

    public static Dictionary<string, (double, double, double, int)> ParseCloneComponents(TextReader fitnessStream)
    {
        var output = new Dictionary<string, (double, double, double, int)>();
        string? firstLine = fitnessStream.ReadLine() ?? throw new Exception("Error in ParseCloneComponents: Fitness file is empty.");
        // Continue past the header
        while (fitnessStream.ReadLine() is { } line)
        {
            var lineSplit = line.Split("\t").Select(s => s.Trim()).ToList();
            string sampleName = lineSplit[0];
            double stressTerm = double.Parse(lineSplit[4], CultureInfo.InvariantCulture.NumberFormat);
            var tsg = double.Parse(lineSplit[5], CultureInfo.InvariantCulture.NumberFormat);
            var og  = double.Parse(lineSplit[6], CultureInfo.InvariantCulture.NumberFormat);
            double tsgogTerm = og + tsg;
            double essTerm = double.Parse(lineSplit[7], CultureInfo.InvariantCulture.NumberFormat);
            // If the file includes data on how many chromosomal events the sample underwent
            int eventCount = lineSplit.Count >= 9 ? int.Parse(lineSplit[8]) : -1;
            output[sampleName] = (stressTerm, tsgogTerm, essTerm, eventCount);
        }
        return output;
    }

    public static Dictionary<string, int> ParseEventCounts(TextReader sampleStream)
    {
        var output = new Dictionary<string, int>();
        string? firstLine = sampleStream.ReadLine() ?? throw new Exception("Error in ParseEventCounts: Sample file is empty.");
        // Continue past the header
        while (sampleStream.ReadLine() is { } line)
        {
            var lineSplit = line.Split("\t").Select(s => s.Trim()).ToList();
            string sampleName = lineSplit[0];
            int eventCount = int.Parse(lineSplit[1]);
            output.Add(sampleName, eventCount);
        }
        return output;
    }

    public static (Dictionary<string, int> chrLengths, Dictionary<string, SexEnum> chrSex) ParseChromosomes(string text)
    {
        IList<string> lines = text.Split("\n");
        Dictionary<string, int> chrLengths = new();
        Dictionary<string, SexEnum> chrSex = new();
        for (var index = 0; index < lines.Count; index++)
        {
            string line = lines[index];
            var lineSplit = line.Split("\t").Select(s => s.Trim()).ToList();
            string chrNo = lineSplit[0];
            int length = int.Parse(lineSplit[1]);
            chrLengths.Add(chrNo, length);
            var sexEnum = GetSexEnum(lines, lineSplit, index);
            chrSex.Add(chrNo, sexEnum);
        }
        return (chrLengths, chrSex);
    }

    public static (Dictionary<string, (int, int)> p, Dictionary<string, (int, int)> q) ParseCentromeres(TextReader centromereFile)
    {   
        centromereFile.ReadLine(); // Skip header
        Dictionary<string, (int, int)> pCentromeres = new();
        Dictionary<string, (int, int)> qCentromeres = new();
        var lastChr = "";
        while (centromereFile.ReadLine() is {} line)
        {
            var lineSplit = line.Split("\t").Select(s => s.Trim()).ToList();
            string chrNo = lineSplit[0];
            int start = int.Parse(lineSplit[1]);
            int end = int.Parse(lineSplit[2]);
            if (chrNo != lastChr)
            {
                lastChr = chrNo;
                pCentromeres.Add(chrNo, (start, end));
            }
            else
            {
                qCentromeres.Add(chrNo, (start, end));
            }
        }
        return (pCentromeres, qCentromeres);
    }

    private static SexEnum GetSexEnum(IList<string> lines, IReadOnlyList<string> lineSplit, int index)
    {
        if (lineSplit.Count > 2)
        {
            string sexString = lineSplit[2];
            return Enum.Parse<SexEnum>(sexString);
        }
        if (index == lines.Count - 1)
        {
            return SexEnum.Male;
        }
        if (index == lines.Count - 2)
        {
            return SexEnum.Female;
        }
        return SexEnum.Both;
    }
    
    public static IEnumerable<StringBuilder> ParseFasta(StreamReader fastaStream)
    {
        StringBuilder? sequence = null;
        while (fastaStream.ReadLine() is { } line)
        {
            if (line.StartsWith(";"))
            {
                continue;
            }
            if (line.StartsWith(">"))
            {
                if (sequence != null)
                {
                    yield return sequence;
                }
                string pattern = @"^>chr([1-9]|1[0-9]|2[0-2]|X|Y)$";
                var match = Regex.Match(line, pattern);
                if (match.Value == "")
                {
                    sequence = null;
                    continue;
                }
                string chrNo = match.Value[1..];
                Console.WriteLine($"Parsing the sequence for chr: " + chrNo);
                // TODO: optimize the string builder
                sequence = new StringBuilder("");
            }
            else if (sequence != null)
            {
                sequence.Append(line);
            }
        }
        if (sequence != null)
        {
            yield return sequence;
        }
    }
}
