// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Collections.Generic;
using System.Linq;

namespace SimChA.DataTypes
{
    public class Chromosome
    {
        private readonly List<Region> _regions;

        public int Length => _regions.Sum(r => r.Length);
        
        public Chromosome(Region initialRegion)
        {
            _regions = new List<Region> { initialRegion };
        }
    }
}