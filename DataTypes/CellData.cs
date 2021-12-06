// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes
{
    public struct CellData
    {
        public int CellId;
        public int ParentId;
        public bool IsAlive;
        public Karyotype Karyotype;
        public float MutationRate;
        public float DivisionRate;

        public CellData(int cellId, int parentId, Karyotype karyotype) : this()
        {
            CellId = cellId;
            ParentId = parentId;
            Karyotype = karyotype;
            IsAlive = true;
            MutationRate = 1f;
            DivisionRate = 1f;
        }
    }
}