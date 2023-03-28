// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record ChrID(ChrNo ChrNo, bool Parent)
{
    private string ParentID() 
        => Parent ?"H1" : "H2";
    
    public override string ToString()
        => $"{ChrNo.ToString()}_{ParentID()}";
}