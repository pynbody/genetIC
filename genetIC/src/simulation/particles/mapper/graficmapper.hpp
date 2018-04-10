#ifndef IC_GRAFICMAPPER_HPP
#define IC_GRAFICMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"

namespace particle {

  namespace mapper {
    using std::endl;
    using std::cerr;


    /*!
         \class GraficMapper
         \brief A particle mapper specifically for Grafic output, where entire grids are written out in order of
         increasing resolution.
        */
    template<typename GridDataType>
    class GraficMapper : public ParticleMapper<GridDataType> {
    public:

      using MapType = ParticleMapper<GridDataType>;
      using typename MapType::T;
      using typename MapType::MapPtrType;
      using typename MapType::GridPtrType;
      using typename MapType::GridType;
      using typename MapType::iterator;
      using typename MapType::ConstGridPtrType;

    protected:
      multilevelcontext::MultiLevelContextInformation<GridDataType> contextInformation;

    public:


      GraficMapper(const multilevelcontext::MultiLevelContextInformation<GridDataType> &context,
                   Coordinate<T> center,
                  size_t extralowres, size_t extrahighres) {
        context.copyContextWithCenteredIntermediate(contextInformation, center, 2, extralowres, extrahighres);
      }

      bool references(GridPtrType grid) const override {
        bool hasReference = false;
        contextInformation.forEachLevel([&grid, &hasReference](const GridType &targetGrid) {
          if (&targetGrid == grid.get())
            hasReference = true;
        });
        return hasReference;
      }

      virtual void debugInfo(std::ostream &s, int level = 0) const override {
        tools::indent(s, level);
        s << "GraficMapper" << endl;
      }

      virtual size_t size() const {
        return contextInformation.getNumCells();
      }

      virtual void flagParticles(const std::vector<size_t> &genericParticleArray) {
        size_t gridStart = 0;
        size_t i = 0;
        contextInformation.forEachLevel([&genericParticleArray, &i, &gridStart](GridType &targetGrid) {
          std::vector<size_t> gridCellArray;
          size_t gridEnd = gridStart + targetGrid.size3;

          // Copy all cells belonging to this grid, offseting the ID appropriately to account for previous grids
          while (i < genericParticleArray.size() && genericParticleArray[i] < gridEnd) {
            gridCellArray.push_back(genericParticleArray[i] - gridStart);
            ++i;
          }

          targetGrid.flagCells(gridCellArray);

          // Update offset
          gridStart = gridEnd;
        });

        if (i != genericParticleArray.size())
          throw std::runtime_error("Ran out of grids when interpreting grafic cell IDs - check IDs?");
      }

      virtual void unflagAllParticles() override {
        contextInformation.forEachLevel([](GridType &targetGrid) {
          targetGrid.unflagAllCells();
        });
      }

      virtual GridPtrType getCoarsestGrid() override {
        return contextInformation.getGridForLevel(0).shared_from_this();
      }

      GridPtrType getFinestGrid() override {
        return contextInformation.getGridForLevel(contextInformation.getNumLevels() - 1).shared_from_this();
      }

      void getFlaggedParticles(std::vector<size_t> &particleArray) const override {
        size_t offset = 0;
        contextInformation.forEachLevel([&particleArray, &offset](const GridType &targetGrid) {
          size_t last = particleArray.size();
          targetGrid.getFlaggedCells(particleArray);

          // Offset according to previous grid sizes
          for (size_t i = last; i < particleArray.size(); ++i) {
            particleArray[i] += offset;
          }

          // Update offset
          offset += targetGrid.size3;
        });
      }

      virtual void
      dereferenceIterator(const iterator * /* *pIterator */, ConstGridPtrType & /*&gp*/,
                          size_t & /*&i*/) const override {
        // Grafic files are written out at the grid level and iterators should not be involved.
        throw std::runtime_error("Iterators are not supported by GraficMapper");
      }

    };
  }
}
#endif
