// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System.Collections.Generic;

namespace SimChA.DataTypes
{
    public class CellList
    {
        private List<CellData> Cells { get; }

        public CellList() : this(new CellData(0, -1, new Karyotype()))
        {
        }
        
        public CellList(CellData initialCell)
        {
            Cells = new List<CellData> { initialCell };
        }

        public void DuplicateCell(int parentId)
        {
            var newCell = new CellData(Cells.Count, parentId, Cells[parentId].Karyotype);
            Cells.Add(newCell);
        }
    }
}