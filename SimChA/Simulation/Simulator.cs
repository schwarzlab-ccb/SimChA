using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
   private readonly CloneList _cloneList;

   public Simulator()
   {
      _cloneList = new CloneList();
   }

   public void Step()
   {
      DivideAndMutate();
      Kill();
   }

   private void Kill()
   {
      throw new NotImplementedException();
   }

   private void DivideAndMutate()
   {
      int cloneCount = _cloneList.Clones.Count;
      for (int i = 0; i < cloneCount; i++)
      {
         var originalClone = _cloneList.Clones[i];
         int newCellsCount = Binomial.Sample(originalClone.DivisionRate, originalClone.AliveCount);
         int newMutantCount = Binomial.Sample(originalClone.MutationRate, newCellsCount); // The existing cells will not mutate
         for (int mutationI = 0; mutationI < newMutantCount; mutationI++)
         {
            var newSubClone = new SubClone(originalClone, _cloneList.Clones.Count);
            var abberation = SelectMutation();
            switch (abberation)
            {
               case AbberationEnum.TailDeletion:
                  newSubClone.Karyotype.ApplyTailDeletion();
                  break;
               case AbberationEnum.Missegregation:
               case AbberationEnum.Duplication:
               case AbberationEnum.Chromothripsis:
               case AbberationEnum.Translocation:
               case AbberationEnum.InternalDeletion:
               case AbberationEnum.Inversion:
               case AbberationEnum.BreakageFusionBridge:
               default:
                  throw new ArgumentOutOfRangeException();
            }
         }
      }
   }

   private AbberationEnum SelectMutation()
   {
      return AbberationEnum.TailDeletion;
   }
}