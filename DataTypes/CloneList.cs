// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes
{
    public class CloneList
    {
        private List<SubClone> Clones { get; }

        public CloneList() : this(new SubClone(0, -1, new Karyotype()))
        {
        }
        
        public CloneList(SubClone initialClone)
        {
            Clones = new List<SubClone> { initialClone };
        }

        public void DuplicateCell(int parentId)
        {
            var subClone = new SubClone(Clones.Count, parentId, Clones[parentId].Karyotype);
            Clones.Add(subClone);
        }
    }
}