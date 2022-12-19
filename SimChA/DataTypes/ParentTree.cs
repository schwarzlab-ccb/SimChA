// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record ParentTree(string RootName, List<TreeNode> Nodes, List<TreeEdge> Edges);