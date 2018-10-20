//
// A simple memmap-enabling class, used in preference to boost because the latter requires building.
//

#ifndef IC_MEMMAP_HPP
#define IC_MEMMAP_HPP

#include <sys/mman.h>
#include <string>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <iostream>

namespace tools {

  class MemMapFileWriter;

  template<typename DataType>
  class MemMapRegion {
  protected:
    char *addr_aligned;
    DataType *addr;
    size_t size_bytes;

    MemMapRegion(int fd, size_t file_offset, size_t n_elements) {
      size_bytes = n_elements*sizeof(DataType);
      ::lseek(fd, file_offset+size_bytes-1, SEEK_SET);
      ::write(fd,"",1);

      size_t npage_offset = file_offset/::getpagesize();
      size_t aligned_offset = npage_offset*::getpagesize();
      size_t byte_page_offset = file_offset-aligned_offset;

      size_bytes+=byte_page_offset;

      addr_aligned = static_cast<char*>(::mmap(nullptr, size_bytes, PROT_READ | PROT_WRITE,
          MAP_SHARED, fd, aligned_offset));

      if(addr_aligned==MAP_FAILED)
        throw std::runtime_error("Failed to create a mem-map for output (reason: "+std::string(::strerror(errno))+")");

      addr = reinterpret_cast<DataType*>(&addr_aligned[byte_page_offset]);

    }

    friend class MemMapFileWriter;

  public:
    ~MemMapRegion() {
      if (addr_aligned != nullptr) {
        msync(addr_aligned, size_bytes, MS_SYNC);
        if(munmap(addr_aligned, size_bytes)!=0) {
          // This probably indicates something has gone catastrophically wrong...
          std::cerr << "ERROR: Failed to delete the mem-map (reason: " << ::strerror(errno) << ")" << std::endl;
          exit(1);
        }
      }
    }

    MemMapRegion(const MemMapRegion & copy) = delete;

    MemMapRegion(MemMapRegion && move) {
      (*this)=std::move(move);
    }

    DataType & operator[](size_t offset) {
      return addr[offset];
    }

    MemMapRegion & operator=(MemMapRegion && move) noexcept {
      this->addr = move.addr;
      this->addr_aligned = move.addr_aligned;
      this->size_bytes = move.size_bytes;
      move.addr_aligned = nullptr;
      return (*this);
    }
  };


  /*!
   \class MemMapFileWriter
   \brief A class that helps to write files with portions of them mem-mapped (particularly helpful for parallel writing)

   It operates on a mixed model where you can sequentially write individual bits of data to the file (using the
   write method), but also get a pointer to the mem-map when convenient (using the getMemMap method).


  */
  class MemMapFileWriter {
  protected:
    int fd;
    size_t offset;

  public:
    MemMapFileWriter() {
      this->fd = -1;
    }

    MemMapFileWriter(const MemMapFileWriter & copy) = delete;

    MemMapFileWriter(MemMapFileWriter && move) {
      (*this)=std::move(move);
    }

    MemMapFileWriter & operator=(MemMapFileWriter &&move) {
      fd = move.fd;
      offset = move.offset;
      move.fd = -1;
      return (*this);
    }



    MemMapFileWriter(std::string filename)  {
      fd = ::open(filename.c_str(), O_RDWR | O_CREAT, (mode_t)0666);
      if(fd==-1)
        throw std::runtime_error("Failed to open file (reason: "+std::string(::strerror(errno))+")");
      ::ftruncate(fd, 0);
      offset = 0;
    }

    ~MemMapFileWriter() {

      if(fd!=-1)
        close(fd);
    }

    //! Write a single item to the file at the current write location
    template<typename DataType>
    void write(const DataType & data) {
      ::write(fd, &data, sizeof(data));
      offset+=sizeof(data);
    }

    //! Write a single item to the file, using Fortran-style size blocks
    template<typename DataType>
    void writeFortran(const DataType &data) {
      constexpr int fortranFieldSize = sizeof(DataType);
      write(fortranFieldSize);
      write(data);
      write(fortranFieldSize);
    }

    //! Get a memory-mapped view of the file at the current write location, with the intention of writing n_elements
    template<typename DataType>
    auto getMemMap(size_t n_elements) {
      auto region = MemMapRegion<DataType>(fd,offset,n_elements);
      offset+=n_elements*sizeof(DataType);
      ::lseek(fd, offset, SEEK_SET);
      return region;
    }

    //! Get a memory-mapped view of the file for writing, and surround it with Fortran-style size blocks
    template<typename DataType>
    auto getMemMapFortran(size_t n_elements) {
      int fortranFieldSize = n_elements*sizeof(DataType);
      write(fortranFieldSize);
      auto region = getMemMap<DataType>(n_elements);
      write(fortranFieldSize);
      return region;
    }

  };
}


#endif //IC_MEMMAP_HPP
