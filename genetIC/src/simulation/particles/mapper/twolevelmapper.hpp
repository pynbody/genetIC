#ifndef IC_TWOLEVELMAPPER_HPP
#define IC_TWOLEVELMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"

namespace particle {

	namespace mapper {
		using std::endl;
		using std::cerr;

		template<typename GridDataType>
		class TwoLevelParticleMapper : public ParticleMapper<GridDataType> {
			/* TwoLevelParticleMapper - the critical class for mapping different grids into an assembled zoom simulation
			 *
			 * This class points to two other mappers, referred to as level 1 and level 2. The level 2 mapper _must_ be
			 * a OneLevelParticleMapper, i.e. it can only point to a single grid. The level 1 mapper on the other hand
			 * can itself be a TwoLevelParticleMapper, allowing multi-level maps to be built through a left-expanding
			 * tree structure.
			 */
		private:

			using MapType = ParticleMapper<GridDataType>;
			using typename MapType::T;
			using typename MapType::MapPtrType;
			using typename MapType::GridPtrType;
			using typename MapType::GridType;
			using typename MapType::iterator;
			using typename MapType::ConstGridPtrType;

			MapPtrType pLevel1;
			MapPtrType pLevel2;

			GridPtrType pGrid1;
			GridPtrType pGrid2;

			size_t n_hr_per_lr;

			size_t totalParticles;
			size_t firstLevel2Particle;

			bool skipLevel1;

			std::vector<size_t> zoomParticleArrayForL1mapper; //< the particles on the level-1 (finest available) grid, which we wish to replace with their zooms
			std::vector<size_t> zoomParticleArrayForL1grid;
			mutable std::vector<size_t> zoomParticleArrayHiresUnsorted; //< the particles on the level-2 grid that are included


			void appendMappedIdsToVector(size_t id0, std::vector<size_t> &ids) const {
				// Finds all the particles at level2 corresponding to the single
				// particle at level1
				T x0, y0, z0;
				std::tie(x0, y0, z0) = pGrid1->getCellCentroid(id0);
				pGrid2->appendIdsInCubeToVector(x0, y0, z0, pGrid1->dx, ids);
			}

			std::vector<size_t> getMappedIds(size_t id0) const {
				std::vector<size_t> ids;
				appendMappedIdsToVector(id0, ids);
				return ids;
			}

			size_t reverseMapId(size_t id2) const {
				auto coord = pGrid2->getCellCentroid(id2);
				return pGrid1->getClosestIdNoWrap(coord);
			}


		public:

			virtual iterator begin(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const override {
				if (zoomParticleArrayForL1mapper.size() == 0)
					return this->end(generator);

				iterator x(this, generator);
				x.extraData.push_back(0); // current position in zoom ID list
				x.extraData.push_back(zoomParticleArrayForL1mapper[0]); // next entry in zoom ID list

				if (!skipLevel1)
					x.subIterators.emplace_back(new iterator(pLevel1->begin(generator)));
				else
					x.subIterators.emplace_back(nullptr);

				x.subIterators.emplace_back(new iterator(pLevel2->begin(generator)));

				syncL2iterator(0, *(x.subIterators.back()));

				return x;
			}


			TwoLevelParticleMapper(MapPtrType &pLevel1,
														 MapPtrType &pLevel2,
														 const std::vector<size_t> &zoomParticlesOnLevel1Grid,
														 bool skipLevel1 = false) :
					pLevel1(pLevel1), pLevel2(pLevel2),
					pGrid1(pLevel1->getFinestGrid()),
					pGrid2(pLevel2->getCoarsestGrid()),
					skipLevel1(skipLevel1),
					zoomParticleArrayForL1grid(zoomParticlesOnLevel1Grid) {
				/* @param pLevel1 - the coarser particle mapper, which can itself have multiple levels
				 * @param pLevel2 - the fine particle mapper, which can point to only one grid (i.e. must be OneLevelParticleMapper)
				 * @param zoomParticlesOnLevel1Grid - the list of particles on the level 1 mapper that will be zoomed. NB these
				 *                                    are specified relative to the finest grid on level 1 (not relative to the
				 *                                    level 1 mapper itself)
				 *
				 *
				 */

				assert(zoomParticleArrayForL1grid.size() > 0);
				assert(typeid(*pLevel2) == typeid(OneLevelParticleMapper<GridDataType>));
				assert(pLevel1->size_gas() == 0);
				assert(pLevel2->size_gas() == 0);

				std::sort(zoomParticleArrayForL1grid.begin(), zoomParticleArrayForL1grid.end());
				pLevel1->unflagAllParticles();
				pGrid1->flagCells(zoomParticleArrayForL1grid);
				pLevel1->getFlaggedParticles(zoomParticleArrayForL1mapper);

				if (zoomParticleArrayForL1grid.size() != zoomParticleArrayForL1mapper.size())
					throw std::runtime_error("The cells to zoom on must all be on the finest level 1 grid");

				n_hr_per_lr = tools::getRatioAndAssertPositiveInteger(pGrid1->dx, pGrid2->dx);
				n_hr_per_lr *= n_hr_per_lr * n_hr_per_lr;

				totalParticles = pLevel1->size() + (n_hr_per_lr - 1) * zoomParticleArrayForL1mapper.size();

				firstLevel2Particle = pLevel1->size() - zoomParticleArrayForL1mapper.size();

				if (skipLevel1) {
					firstLevel2Particle = 0;
					totalParticles = n_hr_per_lr * zoomParticleArrayForL1mapper.size();
				}

				calculateHiresParticleList();

			}

			virtual void debugInfo(std::ostream &s, int level = 0) const override {
				tools::indent(s, level);
				s << "TwoLevelParticleMapper, n_hr_per_lr=" << n_hr_per_lr << ", firstLevel2Particle="
					<< firstLevel2Particle <<
					std::endl;
				tools::indent(s, level);
				s << "                      , zoom.size=" << zoomParticleArrayForL1mapper.size() << ", zoomed.size=" <<
					zoomParticleArrayHiresUnsorted.size() << std::endl;
				if (skipLevel1) {
					tools::indent(s, level);
					s << "low-res part will be skipped but notionally is:" << endl;
				}
				pLevel1->debugInfo(s, level + 1);
				tools::indent(s, level);
				s << "TwoLevelParticleMapper continues with high-res particles:" << std::endl;
				pLevel2->debugInfo(s, level + 1);
				tools::indent(s, level);
				s << "TwoLevelParticleMapper ends" << std::endl;
			}

			bool references(GridPtrType grid) const override {
				return pLevel1->references(grid) || pLevel2->references(grid);
			}

			virtual void unflagAllParticles() override {
				pLevel1->unflagAllParticles();
				pLevel2->unflagAllParticles();
			}

			void flagParticles(const std::vector<size_t> &orderedParticleIndices) override {

				std::vector<size_t> level1particles;
				std::vector<size_t> level2particles;

				if (pLevel1->size() < zoomParticleArrayForL1mapper.size())
					throw std::runtime_error("Zoom particle list is longer than the grid it refers to");

				size_t i0 = pLevel1->size() - zoomParticleArrayForL1mapper.size();
				size_t lr_index;
				size_t lr_mapper_particle;

				size_t last_val = 0;
				size_t numSkippedParticlesLR = 0;
				size_t nZoomParticles = zoomParticleArrayForL1mapper.size();

				for (size_t i : orderedParticleIndices) {
					if (i < last_val)
						throw std::runtime_error("The particle list must be in ascending order");
					last_val = (i);

					if (i < i0) {
						// particle in low res region. We need to count how many skipped
						// particles there are to work out the address in the original
						// file ordering. NB this approach assumes the input particle array is
						// in ascending order (and the zoomParticleArray).
						while (numSkippedParticlesLR < nZoomParticles &&
									 zoomParticleArrayForL1mapper[numSkippedParticlesLR] <= i + numSkippedParticlesLR)
							++numSkippedParticlesLR;

						level1particles.push_back(i + numSkippedParticlesLR);

					} else {
						// particle in high res region
						lr_index = (i - i0) / n_hr_per_lr;

						if (lr_index >= zoomParticleArrayForL1mapper.size())
							throw std::runtime_error("Particle ID out of range");

						// find the low-res particle. Note that we push it onto the list
						// without searching to see if it's already there, i.e. we assume that
						// duplicates don't matter.

						lr_mapper_particle = zoomParticleArrayForL1mapper[lr_index];
						level1particles.push_back(lr_mapper_particle);

						// get all the HR particles
						std::vector<size_t> hr_particles = getMappedIds(zoomParticleArrayForL1grid[lr_index]);
						assert(hr_particles.size() == n_hr_per_lr);

						// work out which of these this particle must be and push it onto
						// the HR list
						size_t offset = (i - i0) % n_hr_per_lr;
						level2particles.push_back(hr_particles[offset]);
					}

				}

				// Some routines need these lists to be sorted (e.g. getFlaggedParticles below,
				// or in principle there could be an underlying further processing step.
				// So it's time to do it...
				std::sort(level1particles.begin(), level1particles.end());
				std::sort(level2particles.begin(), level2particles.end());

				// Get the grids to interpret their own particles. This almost
				// certainly just involves making a copy of our list (notwithstanding
				// point made above about a possible underlying processing step). In fact, we're
				// done with our list so they might as well steal the data.
				pLevel1->flagParticles(std::move(level1particles));
				pLevel2->flagParticles(std::move(level2particles));

			}


			void getFlaggedParticles(std::vector<size_t> &particleArray) const override {

				// translate level 1 particles - just need to exclude the zoomed particles
				std::vector<size_t> grid1particles;
				pLevel1->getFlaggedParticles(grid1particles);

				size_t zoom_i = 0;
				size_t len_zoom = zoomParticleArrayForL1mapper.size();
				for (const size_t &i_lr : grid1particles) {
					while (zoom_i < len_zoom && zoomParticleArrayForL1mapper[zoom_i] < i_lr)
						++zoom_i;

					if (zoom_i == len_zoom || zoomParticleArrayForL1mapper[zoom_i] != i_lr) {
						// not a zoom particle: record it in the low res region
						particleArray.push_back(i_lr - zoom_i);
					}
				}

				std::vector<size_t> grid2particles;
				pLevel2->getFlaggedParticles(grid2particles);

				// get the ordering of the level 2 particles
				std::vector<size_t> sortIndex = tools::argsort(zoomParticleArrayHiresUnsorted);


				size_t zoomed_i = 0;
				size_t len_zoomed = zoomParticleArrayHiresUnsorted.size();
				for (const size_t &i_hr : grid2particles) {
					// find the i_hr in the zoomed particle list
					while (zoomed_i < len_zoomed && zoomParticleArrayHiresUnsorted[sortIndex[zoomed_i]] < i_hr)
						++zoomed_i;

					// If the marked particle is not actually in the output list, ignore it.
					//
					// Older versions of the code throw an exception instead
					if (zoomed_i == len_zoomed || zoomParticleArrayHiresUnsorted[sortIndex[zoomed_i]] != i_hr)
						continue;

					particleArray.push_back(sortIndex[zoomed_i] + firstLevel2Particle);
				}

				std::sort(particleArray.begin(), particleArray.end());

			}

			virtual GridPtrType getCoarsestGrid() override {
				return pLevel1->getCoarsestGrid();
			}

			GridPtrType getFinestGrid() override {
				return pLevel2->getFinestGrid();
			}

			virtual size_t size() const override {
				return totalParticles;
			}

			std::pair<MapPtrType, MapPtrType> addGas(T massRatio, const std::vector<GridPtrType> &toGrids) override {
				bool newskip = skipLevel1;

				auto newLevel1 = pLevel1->addGas(massRatio, toGrids);
				auto newLevel2 = pLevel2->addGas(massRatio, toGrids);

				auto gasSubLevel1 = newLevel1.first;
				auto gasSubLevel2 = newLevel2.first;
				auto dmSubLevel1 = newLevel1.second;
				auto dmSubLevel2 = newLevel2.second;

				decltype(gasSubLevel1) newGasMap;
				decltype(gasSubLevel1) newDmMap;

				if (gasSubLevel1 == nullptr) {
					gasSubLevel1 = pLevel1;
					newskip = true;
				}

				if (gasSubLevel2 != nullptr)
					newGasMap = std::make_shared<TwoLevelParticleMapper<GridDataType>>(
							gasSubLevel1, gasSubLevel2, zoomParticleArrayForL1mapper,
							newskip);
				else
					newGasMap = nullptr;

				newDmMap = std::make_shared<TwoLevelParticleMapper<GridDataType>>(
						dmSubLevel1, dmSubLevel2, zoomParticleArrayForL1mapper,
						skipLevel1);

				return std::make_pair(newGasMap, newDmMap);


			}

			MapPtrType superOrSubSampleDM(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {

				auto ssub1 = pLevel1->superOrSubSampleDM(ratio, toGrids, super);
				auto ssub2 = pLevel2->superOrSubSampleDM(ratio, toGrids, super);

				// Work out the new list of particles to zoom on
				decltype(zoomParticleArrayForL1mapper) newZoomParticles;
				ssub1->unflagAllParticles();
				pLevel1->flagParticles(zoomParticleArrayForL1mapper);
				ssub1->getFlaggedParticles(newZoomParticles);

				return std::make_shared<TwoLevelParticleMapper<GridDataType>>(
						ssub1, ssub2, newZoomParticles,
						skipLevel1);
			}


		protected:

			void syncL2iterator(const size_t &next_zoom, iterator &level2iterator) const {

				if (zoomParticleArrayHiresUnsorted[next_zoom] > level2iterator.i) {
					level2iterator += zoomParticleArrayHiresUnsorted[next_zoom] - level2iterator.i;
				} else {
					level2iterator -= level2iterator.i - zoomParticleArrayHiresUnsorted[next_zoom];
				}

				assert(level2iterator.i == zoomParticleArrayHiresUnsorted[next_zoom]);
			}

			void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const override {
				tools::indent(s, n);
				s << "i=" << pIterator->i << " into iterator for TwoLevelParticleMapper " << endl;
				iterator &level1iterator = *(pIterator->subIterators[0]);
				iterator &level2iterator = *(pIterator->subIterators[1]);
				const size_t &i = pIterator->i;
				const size_t &next_zoom = pIterator->extraData[0];
				const size_t &next_zoom_index = pIterator->extraData[1];

				if (i >= firstLevel2Particle) {
					tools::indent(s, n);
					s << "Inside level 2 particles" << endl;
					level2iterator.debugInfo(s, n + 1);
				} else {
					tools::indent(s, n);
					s << "Inside level 1 particles" << endl;
					tools::indent(s, n);
					s << "next_zoom = " << next_zoom << "; next_zoom_index = " << next_zoom_index << endl;
					level1iterator.debugInfo(s, n + 1);
				}
			}

			virtual void incrementIterator(iterator *pIterator) const override {

				// set up some helpful shortcut references
				auto &extraData = pIterator->extraData;
				iterator &level1iterator = *(pIterator->subIterators[0]);
				iterator &level2iterator = *(pIterator->subIterators[1]);
				size_t &i = pIterator->i;

				size_t &next_zoom = extraData[0];
				size_t &next_zoom_index = extraData[1];

				// increment the boring-old-counter!
				++i;

				// now work out what it actually points to...

				if (i >= firstLevel2Particle) {
					// do zoom particles

					if (i == firstLevel2Particle)
						next_zoom = 0;
					else
						++next_zoom;

					syncL2iterator(next_zoom, level2iterator);


				} else {
					// do normal particles
					//
					// N.B. should never reach here if level1iterator is NULL

					assert((&level1iterator) != nullptr);

					++level1iterator;

					while (level1iterator.i == next_zoom_index) {
						// we DON'T want to return this particle.... it will be 'zoomed'
						// later on... ignore it
						++level1iterator;

						// load the next zoom particle
						++next_zoom;
						if (next_zoom < zoomParticleArrayForL1mapper.size())
							next_zoom_index = zoomParticleArrayForL1mapper[next_zoom];
						else
							next_zoom_index = size() + 1; // i.e. there isn't a next zoom index!

						// N.B. now it's possible we our 'next' particle is also to be zoomed on,
						// so the while loop goes back round to account for this
					}


				}

			}


			virtual void incrementIteratorBy(iterator *pIterator, size_t increment) const {
				// could be optimized:
				for (size_t i = 0; i < increment; ++i)
					incrementIterator(pIterator);

			}

			virtual void
			dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
				if (pIterator->i >= firstLevel2Particle)
					pIterator->subIterators[1]->deReference(gp, i);
				else
					pIterator->subIterators[0]->deReference(gp, i);
			}

			void calculateHiresParticleList() const {

				zoomParticleArrayHiresUnsorted.reserve(zoomParticleArrayForL1mapper.size() * n_hr_per_lr);

				for (size_t i : zoomParticleArrayForL1grid) {
					appendMappedIdsToVector(i, zoomParticleArrayHiresUnsorted);

				}

				if (!pLevel2->supportsReverseIterator()) {
					// underlying map can't cope with particles being out of order - sort them
					std::sort(zoomParticleArrayHiresUnsorted.begin(), zoomParticleArrayHiresUnsorted.end());
				}

			}


		};
	}
}
#endif
