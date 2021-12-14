using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

namespace SimChA.Simulation;

public class Simulator
{
   private CloneList _cloneList;
   private Random _random;

   public Simulator()
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
         ChrMutations.DeleteRegion(randomChromosome, 0, deletionLength);
      }
      else
      {
         ChrMutations.DeleteRegion(randomChromosome, chromLength - deletionLength, chromLength);
      }
   }
}