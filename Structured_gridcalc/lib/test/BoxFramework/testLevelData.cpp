#include <iostream>
#include <iomanip>

#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"

int main(const int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;

//--Tests

  Box domain(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(7, 7, 7)));

#if 1
  // Using box size 4 should result in sizes (4, 4)
  DisjointBoxLayout dbl(domain, 4*IntVect::Unit);
  if (dbl.size() != ((2*IntVect::Unit).product())) ++status;
  // First box should have maxsize
  if (dbl.getLinear(0).box.size() != ((4*IntVect::Unit).product())) ++status;
  // Second box should have maxsize
  if (dbl.getLinear(1).box.size() != ((4*IntVect::Unit).product())) ++status;
  // Last box should have maxsize
  if (dbl.getLinear(dbl.size()-1).box.size() !=
      ((4*IntVect::Unit).product())) ++status;

  const int numBox = D_TERM(2, *2, *2);

  // Test construction and define of simple level data of type LDTest
  {
    class LDTest
    {
    public:
      using value_type = int;
      LDTest()
        :
        m_boxSize(0),
        m_ncomp(0)
        { }
      void define(const Box& a_box, const int a_ncomp)
        {
          m_boxSize = a_box.size();
          m_ncomp = a_ncomp;
        }
      int m_boxSize;
      int m_ncomp;
    };
    LevelData<LDTest> lvldata(dbl, 2, 0);
    for (int c = 0; c != numBox; ++c)
      {
        const LDTest& ldt = lvldata.getLinear(c);
        if (ldt.m_boxSize != (4*IntVect::Unit).product()) ++status;
        if (ldt.m_ncomp != 2) ++status;
      }
    
    lvldata.define(dbl, 1, 1);
    for (int c = 0; c != numBox; ++c)
      {
        const LDTest& ldt = lvldata.getLinear(c);
        if (ldt.m_boxSize != (6*IntVect::Unit).product()) ++status;
        if (ldt.m_ncomp != 1) ++status;
      }
  
    // Test DataIterator
    for (DataIterator dit(dbl); dit.ok(); ++dit)
      {
        const LDTest& ldt = lvldata[dit];
        if (ldt.m_boxSize != (6*IntVect::Unit).product()) ++status;
        if (ldt.m_ncomp != 1) ++status;
      }
  }
  
  // Test settting of all values
  LevelData<BaseFab<Real> > lvldata(dbl, 2, 1);
  lvldata.setVal(2.5);

  for (int c = 0; c != numBox; ++c)
  {
    const BaseFab<Real>& fab = lvldata.getLinear(c);
    for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
      {
        if (fab(*bit, 0) != 2.5 || fab(*bit, 1) != 2.5) ++status;
      }
  }

  // Test settting of single component
  lvldata.setVal(1, -2.5);

  for (int c = 0; c != numBox; ++c)
    {
      const BaseFab<Real>& fab = lvldata.getLinear(c);
      for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
        {
          if (fab(*bit, 0) != 2.5 || fab(*bit, 1) != -2.5) ++status;
        }
    }

  // Test setting values per BaseFab
  int c = 1;
  for (DataIterator dit(dbl); dit.ok(); ++dit, ++c)
    {
      lvldata[dit].setVal(0, (Real)c);
      lvldata[dit].setVal(1, (Real)(-c));
    }

  c = 0;
  D_INVTERM(for (int i = 0; i != 2; ++i),
            for (int j = 0; j != 2; ++j),
            for (int k = 0; k != 2; ++k))
    {
      const int boxSize = (4*IntVect::Unit).product();
      const Box& box = dbl.getLinear(c).box;
      if (box.size() != boxSize) ++status;
      const BaseFab<Real>& fab = lvldata.getLinear(c);
      ++c;
      for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
        {
          if (fab(*bit, 0) != ((Real)c) || fab(*bit, 1) != ((Real)(-c)))
            ++status;
        }
    }
#endif

#if 1
  // Test exchange
  if (verbose) std::cout << "Testing exchange\n";
  Copier copier;
  copier.defineExchangeLD<BaseFab<Real> >(lvldata);
  // copier.defineExchangeDBL<Real>(dbl, 1, 0, 2);
  lvldata.exchange(copier);
  {
    c = 0;
    Box domainIVBox(IntVect::Zero, IntVect::Unit);  // Two boxes in each Dir.
    for (BoxIterator domainIVit(domainIVBox); domainIVit.ok();
         ++domainIVit, ++c)
      {
        const BaseFab<Real>& fab = lvldata.getLinear(c);
        if (verbose) std::cout << "c: " << c << std::endl;
        const Box& box = dbl.getLinear(c).box;
        if (verbose) std::cout << "The box: " << box << std::endl;
        // The IntVect from the iterator is the direction to the neighbour
        for (BoxIterator nbrIVit(Box(-IntVect::Unit, IntVect::Unit));
             nbrIVit.ok(); ++nbrIVit)
          {
            // IntVect describing the location of the neighbour in domainIVBox
            const IntVect nbrIV = *domainIVit + *nbrIVit;
            if ((*nbrIVit != IntVect::Zero) && domainIVBox.contains(nbrIV))
              {
                int cNbr = D_TERM(nbrIV[0]*1, + nbrIV[1]*2, + nbrIV[2]*4);
                const Box& nbrBox = dbl.getLinear(cNbr).box;
                Box region(box);
                region.shift(*nbrIVit);
                region &= nbrBox;
                if (verbose)
                  {
                    std::cout << "  Neighbour dir: " << *nbrIVit << std::endl;
                    std::cout << "  Region: " << region << std::endl;
                  }
                // The value to expect in this region
                const int val = 1 + cNbr;
                for (BoxIterator bit(region); bit.ok(); ++bit)
                  {
                    if (fab(*bit, 0) != ((Real)val) ||
                        fab(*bit, 1) != ((Real)(-val))) ++status;
                  }
              }
          }
      }
  }
#endif

//--Output status

  if (verbose)
    {
      std::cout << "Status: " << status << std::endl;
    }
  const char* const testName = "testLevelData";
  const char* const statLbl[] = {
    "failed",
    "passed"
  };
  std::cout << std::left << std::setw(40) << testName
            << statLbl[(status == 0)] << std::endl;
  return status;
}
