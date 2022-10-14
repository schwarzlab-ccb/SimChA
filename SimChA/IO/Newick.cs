using SimChA.DataTypes;
using SimChA.Simulation;
using System.Text.RegularExpressions;

namespace SimChA.IO;

public static class Newick
{
    private static Clone CreateNodes(string newickNode, int parentId, bool isFemale)
    {
        // TODO: all nodes should be considered alive / don't parse population size
        string[] cloneString = newickNode.Split(':');
        // TODO: split below in individual assignments
        var clone = new Clone(int.Parse(cloneString[0].Split('-')[0]), parentId, "1", int.Parse(cloneString[1]), new Karyotype(isFemale));
        return clone;
    }


    public static List<Clone> ParseNewick(string newickString, bool isFemale)
    {
        List<Clone> clones = new();
        string regexPattern = @"(?<closeChildren>[(])|" +
                              @"(?<openChildren>[)])|" +
                              @"(?<branchLength>:+[0-9a-zA-Z-_.]+)|" +
                              @"(?<nodeName>[0-9a-zA-Z-_.]+)|" +
                              @"(?<nextNode>[,])|" +
                              @"(?<root>[;])";
                              /*@"(\()|(\))|([0-9a-zA-Z-_.]+)|(\:)|(\\n)|(\;)|(\,)"*/ 
        RegexOptions regexOptions = RegexOptions.RightToLeft | RegexOptions.IgnorePatternWhitespace;
        if(newickString == "")
        {
            return clones;
        }
        var matches = Regex.Matches(newickString, regexPattern, regexOptions);
        if(!matches.Any() || matches[0].Value != ";") 
        {
            throw new Exception("Newick file is not in the right format");
        }
        var branchLength = checkBranchLength(matches);
        //Create Clones, reversed order to start with root
        var parentIds = new List<int> {-1};
        foreach(Match match in matches)
        {
            
            switch(match.Value)
            {
                case ";":
                    //create root node
                    clones.Add(createClone(clones.Count(), parentIds.Last(), match.NextMatch(), 
                        match.NextMatch().NextMatch(), isFemale, branchLength));
                    parentIds.Add(clones[clones.Count()-1].CloneId);
                    break;
                case ")":
                    //create new children
                    parentIds.Add(clones[clones.Count()-1].CloneId);
                    clones.Add(createClone(clones.Count(), parentIds.Last(), match.NextMatch(), 
                        match.NextMatch().NextMatch(), isFemale, branchLength));
                    clones[parentIds.Last()].ChildrenIDs.Add(clones.Count()-1);
                    break;
                case "(":
                    //remove parent from parentList and add childrenID to parent
                    parentIds.RemoveAt(parentIds.Count()-1);
                    break;
                case ",":
                    //create new child
                    clones.Add(createClone(clones.Count(), parentIds.Last(), match.NextMatch(), 
                        match.NextMatch().NextMatch(), isFemale, branchLength));
                    clones[parentIds.Last()].ChildrenIDs.Add(clones.Count()-1);
                    break;
            }
        }
        if(!clones.Any())
        {
            throw new Exception("No clones found in newick file. Might not be the right format.");
        }
        return clones;
    }
    //create clone from newick Match
    private static Clone createClone(int id, int parentId, Match branchLengthMatch, Match nameMatch, bool isFemale, bool branchLength)
    {
        string nameClone = nameMatch.Groups["nodeName"].Value != "" ? nameMatch.Value : 
            branchLengthMatch.Groups["nodeName"].Value != "" ? branchLengthMatch.Value : id.ToString();
        int mutCount = branchLengthMatch.Groups["branchLength"].Value != "" ? 
            (int)Math.Ceiling(float.Parse(branchLengthMatch.Value.Remove(0,1))) : branchLength ? 0 : 1;
        var clone = new Clone(id, parentId, nameClone, mutCount, new Karyotype(isFemale));
        return clone;
    }
    //check for branch-length in newick file
    private static bool checkBranchLength(MatchCollection matches)
    {
        bool branchLength = false;
        foreach (Match match in matches)
        {
            if(match.Groups["branchLength"].Value != "")
            {
                branchLength = true;
                break;
            }
        }
        Console.Write(branchLength ? "" : "No branch-lengths were found, using 1 as branch-length.");
        return branchLength;
    }
    
    public static List<Clone> ParseNewickString(string[] newickString, bool isFemale)
    {
        List<Clone> clones = new();
        var parentIds = new List<int> { -1 };
        bool rootSet = false;
        // TODO: Multiple code repetitions below, fix
        for (int i = 0; i < newickString.Length; i++)
        {
            switch (newickString[i])
            {
                case "(":
                    if (newickString[i - 1] == "")
                    {
                        parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                        break;
                    }

                    clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                    parentIds = parentIds.Where(p => p != parentIds.Last()).ToList();
                    break;
                case ")":
                    if (rootSet)
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                    }

                    break;
                case ",":
                    if (!rootSet)
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                        parentIds.Add(int.Parse(newickString[i - 1].Split('-')[0]));
                        rootSet = true;
                    }
                    else
                    {
                        clones.Add(CreateNodes(newickString[i - 1], parentIds.Last(), isFemale));
                    }

                    break;
            }
        }
        return clones;
    }
}