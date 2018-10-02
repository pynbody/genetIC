//
// A simple memmap-enabling class, used in preference to boost because the latter requires building.
//

#ifndef IC_MEMMAP_HPP
#define IC_MEMMAP_HPP

#include <sys/mman.h>
#include <string>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

namespace tools {


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
    size_t file_length;
    char *addr;

  public:
    MemMapFileWriter(const MemMapFileWriter & copy) = delete;

    MemMapFileWriter(MemMapFileWriter && move) {
      this->fd = move.fd;
      this->addr = move.addr;
      this->offset = move.offset;
      this->file_length = move.file_length;
      move.fd = -1;
      move.addr = nullptr;
    }

    MemMapFileWriter(std::string filename, size_t file_length) : file_length(file_length) {
      fd = open(filename.c_str(), O_RDWR | O_CREAT, (mode_t)0666);
      if(fd==-1)
        throw std::runtime_error("Failed to open file (reason "+std::to_string(errno)+")");
      lseek(fd, file_length-1, SEEK_SET);
      ::write(fd,"",1);
      addr = static_cast<char*>(mmap(nullptr, file_length, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
      if(addr==MAP_FAILED)
        throw std::runtime_error("Failed to create a mem-map for output (reason "+std::to_string(errno)+")");

      offset = 0;
    }

    ~MemMapFileWriter() {
      if(addr!=nullptr) {
        msync(static_cast<void *>(addr), file_length, MS_SYNC);
        munmap(static_cast<void *>(addr), file_length);
      }
      if(fd!=-1)
        close(fd);
    }

    //! Write a single item to the file at the current write location
    template<typename DataType>
    void write(const DataType & data) {
      assert(offset+sizeof(DataType)<=file_length);
      *(reinterpret_cast<DataType*>(&(addr[offset]))) = data;
      offset+=sizeof(data);
    }

    //! Get a memory-mapped view of the file at the current write location, with the intention of writing n_elements
    template<typename DataType>
    DataType *getMemMap(size_t n_elements) {
      assert(offset+n_elements*sizeof(DataType)<=file_length);
      DataType* region=reinterpret_cast<DataType*>(&(addr[offset]));
      offset+=n_elements*sizeof(DataType);
      return region;
    }

  };
}


#endif //IC_MEMMAP_HPP
