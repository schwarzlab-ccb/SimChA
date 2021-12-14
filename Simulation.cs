using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA;

public class Simulation
{
   private CloneList _cloneList;
   private Random _random;

   public Simulation()
   {
      _cloneList = new CloneList();
      _random = new Random();
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
                  ApplyTailDeletion(newSubClone.Karyotype);
                  break;
               default:
                  break;
            }
         }
      }
   }

   private AbberationEnum SelectMutation()
   {
      return AbberationEnum.TailDeletion;
   }

   private void ApplyTailDeletion(Karyotype karyotype)
   {
      var randomChromosome = karyotype.Chromosomes.Shuffle().First();
      int chromLength = randomChromosome.Length;
      int deletionLength = _random.Next(0, chromLength);
      if (_random.CoinFlip())
      {
         randomChromosome.DeleteRegion(0, deletionLength);
      }
      else
      {
         randomChromosome.DeleteRegion(chromLength - deletionLength, chromLength);
      }
   }
}