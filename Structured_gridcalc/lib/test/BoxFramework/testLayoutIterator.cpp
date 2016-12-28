#include <iostream>
#include <iomanip>

#include "LayoutIterator.H"

int main(const int argc, const char* argv[])
{
  const bool verbose = ((argc == 2) && (std::strcmp(argv[1], "-v") == 0));
  int status = 0;

//--Tests

  // To make life a little more interesting, use a more complicated domain
  // The size is still 12 in each direction
  IntVect domainLo(D_DECL(1, -1, 1));
  Box domain(domainLo, domainLo + 11*IntVect::Unit);

//#if 1
  // Using box size 4
  DisjointBoxLayout dbl(domain, 4*IntVect::Unit);
  if (dbl.size() != ((3*IntVect::Unit).product())) ++status;
  // First box should have maxsize
  if (dbl.getLinear(0).box.size() != ((4*IntVect::Unit).product())) ++status;
  // Last box should have maxsize
  if (dbl.getLinear(dbl.size()-1).box.size() != ((4*IntVect::Unit).product()))
    ++status;

  // If the boxes in the domain are represented by IntVects, the layoutIterator
  // should follow the boxIterator
  {
    // Test LayoutIterator
    BoxIterator boxit(Box(IntVect::Zero, 2*IntVect::Unit));
    for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++boxit)
      {
        CH_assert(boxit.ok());
        const IntVect& ivBox = *boxit;
        const IntVect lo = 4*ivBox;
        IntVect hi = (4*(ivBox + IntVect::Unit) - IntVect::Unit);
        Box testBox(lo + domain.loVect(), hi + domain.loVect());
        if (verbose)
          {
            std::cout << "testBox: " << testBox << std::endl;
            std::cout << "Lit    : " << dbl[lit] << std::endl;
          }
        if (dbl[lit] != testBox) ++status;
      }
  
  //#if 1
    // Test DataIterator which is same as LayoutIterator for single processor
    boxit = BoxIterator(Box(IntVect::Zero, 2*IntVect::Unit));  // Reset
    for (DataIterator dit(dbl); dit.ok(); ++dit, ++boxit)
      {
        CH_assert(boxit.ok());
        const IntVect& ivBox = *boxit;
        const IntVect lo = 4*ivBox;
        IntVect hi = (4*(ivBox + IntVect::Unit) - IntVect::Unit);
        Box testBox(lo + domain.loVect(), hi + domain.loVect());
        if (verbose)
          {
            std::cout << "testBox: " << testBox << std::endl;
            std::cout << "Dit    : " << dbl[dit] << std::endl;
          }
        if (dbl[dit] != testBox) ++status;
      }
  }

//--Test the neighbor iterator

  {
    Box domainIVBox(IntVect::Zero, 2*IntVect::Unit);
    BoxIterator domainIVit(domainIVBox);
    for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++domainIVit)
      {
        const IntVect& ivBase = *domainIVit;  // IV nbrs are centered around
        Box nbrBox(ivBase, ivBase);
        nbrBox.grow(1);
        nbrBox &= domainIVBox;                // Box describing neighbours
        BoxIterator boxit(nbrBox);
        if (*boxit == ivBase) ++boxit;
        for (NeighborIterator nbrit(lit); nbrit.ok(); ++nbrit)
          {
            const IntVect& ivNbrBox = *boxit;
            const IntVect lo = 4*ivNbrBox;
            IntVect hi = (4*(ivNbrBox + IntVect::Unit) - IntVect::Unit);
            Box testNbrBox(lo + domain.loVect(), hi + domain.loVect());
            if (dbl[nbrit] != testNbrBox) ++status;
            if (verbose && status)
              {
                std::cout << "This : " << dbl[lit] << std::endl;
                std::cout << "Nbr  : " << dbl[nbrit] << std::endl;
                std::cout << "test : " << testNbrBox << std::endl;
                std::cout << "NbrIV: " << ivNbrBox << std::endl;
              }
            ++boxit;
            if (*boxit == ivBase) ++boxit;
          }
      }
  }


//--Test the neighbor iterator with trimming faces, edges, and corners
//--separately

  for (int iTrim = 1; iTrim != 4; ++iTrim)
  {
    const unsigned trim = (1<<iTrim);
    Box domainIVBox(IntVect::Zero, 2*IntVect::Unit);
    BoxIterator domainIVit(domainIVBox);
    for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++domainIVit)
      {
        const IntVect& ivBase = *domainIVit;  // IV nbrs are centered around
        Box nbrBox(ivBase, ivBase);
        nbrBox.grow(1);
        nbrBox &= domainIVBox;                // Box describing neighbours
        BoxIterator boxit(nbrBox);
        while (boxit.ok() &&
               (1 << (*boxit - ivBase).norm1()) & (1 + trim)) ++boxit;
        for (NeighborIterator nbrit(lit, trim); nbrit.ok(); ++nbrit)
          {
            const IntVect& ivNbrBox = *boxit;
            const IntVect lo = 4*ivNbrBox;
            IntVect hi = (4*(ivNbrBox + IntVect::Unit) - IntVect::Unit);
            Box testNbrBox(lo + domain.loVect(), hi + domain.loVect());
            if (dbl[nbrit] != testNbrBox) ++status;
            if (verbose && status)
              {
                std::cout << "This : " << dbl[lit] << std::endl;
                std::cout << "Nbr  : " << dbl[nbrit] << std::endl;
                std::cout << "test : " << testNbrBox << std::endl;
                std::cout << "NbrIV: " << ivNbrBox << std::endl;
              }
            ++boxit;
            while (boxit.ok() && (1 << (*boxit - ivBase).norm1()) & (1 + trim))
              ++boxit;
          }
      }
  }
    #if 1
//--Test the periodic iterator

  {
    Box domainIVBox(IntVect::Zero, 2*IntVect::Unit);
    Box perDomain(domainIVBox);
    perDomain.grow(1);
    Box perDomainSide[g_SpaceDim][2];
    for (int dir = 0; dir != g_SpaceDim; ++dir)
      {
        perDomainSide[dir][0] = perDomain;
        perDomainSide[dir][0].grow(-1, dir);
        perDomainSide[dir][1] = perDomainSide[dir][0];
        perDomainSide[dir][0].adjBox(1, dir, -1);
        perDomainSide[dir][1].adjBox(1, dir, 1);
      }
    BoxIterator domainIVit(domainIVBox);
    for (LayoutIterator lit(dbl); lit.ok(); ++lit, ++domainIVit)
      {
        const IntVect& ivBase = *domainIVit;  // IV nbrs are centered around
        Box nbrBox(ivBase, ivBase);
        nbrBox.grow(1);
        BoxIterator boxit(nbrBox);
        if (*boxit == ivBase) ++boxit;
        while (boxit.ok() && domainIVBox.contains(*boxit)) ++boxit;
        
        for (PeriodicIterator perit(lit, 0, 7); perit.ok(); ++perit)
          {
            IntVect ivNbrBox = *boxit;
            for (int dir = 0; dir != g_SpaceDim; ++dir)
              {
                if (perDomainSide[dir][0].contains(ivNbrBox))
                  ivNbrBox[dir] += 3;
                if (perDomainSide[dir][1].contains(ivNbrBox))
                  ivNbrBox[dir] -= 3;
              }
            const IntVect lo = 4*ivNbrBox;
            IntVect hi = (4*(ivNbrBox + IntVect::Unit) - IntVect::Unit);
            Box testPerBox(lo + domain.loVect(), hi + domain.loVect());
            if(dbl[perit] != testPerBox) ++status;
            if (verbose && status)
              {
                std::cout << "This     : " << dbl[lit] << std::endl;
                std::cout << "Periodic : " << dbl[perit] << std::endl;
                std::cout << "test     : " << testPerBox << std::endl;
                std::cout << "NbrIV    : " << ivNbrBox << std::endl;
              }
            ++boxit;
            while (boxit.ok() && domainIVBox.contains(*boxit)) ++boxit;
          }
      }
  }
#endif

//--Output status

  if (verbose)
    {
      std::cout << "Status: " << status << std::endl;
    }
  const char* const testName = "testLayoutIterator";
  const char* const statLbl[] = {
    "failed",
    "passed"
  };
  std::cout << std::left << std::setw(40) << testName
            << statLbl[(status == 0)] << std::endl;
  return status;
}
