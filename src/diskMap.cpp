#include <maki/diskMap.hpp>
#include <maki/utils.hpp>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <vector>

namespace pmap
{
    
    namespace alloc
    {

        TbbBlockAllocator::TbbBlockAllocator(size_t blockSize) : _bsize(blockSize) {}

        memory::memory_block TbbBlockAllocator::allocate_block()
        {
            char *mem = oneapi::tbb::tbb_allocator<char>{}.allocate(_bsize);
            return { reinterpret_cast<void *>(mem), _bsize };
        }

        void TbbBlockAllocator::deallocate_block(memory::memory_block b)
        {
            oneapi::tbb::tbb_allocator<char>{}.deallocate(
                reinterpret_cast<char *>(b.memory),
                b.size);
        }

        size_t TbbBlockAllocator::next_block_size() const noexcept
        { return _bsize; }

        MMapBlockAllocator::File::File(std::string &&file_, size_t bytes_)
            : file(initialiseFile(std::move(file_), bytes_))
            , mmap(file)
            , bytes(bytes_)
        {
        }

        std::string &&MMapBlockAllocator::File::initialiseFile(std::string &&file_, size_t bytes_)
        {
            std::ofstream{file_};
            std::filesystem::resize_file(file_, bytes_);
            return std::move(file_);
        }

        void MMapBlockAllocator::File::sync()
        {
            std::error_code error;
            mmap.sync(error);
            if (error)
                throw std::runtime_error(
                    "Failed to sync mmap to disk (" + std::to_string(error.value()) + "): " + error.message());
        }

        MMapBlockAllocator::File::~File()
        {
            sync();
            mmap.unmap();
            std::filesystem::remove(file);
        }

        MMapBlockAllocator::MMapBlockAllocator(size_t blockSize, std::string_view base)
            : _bsize(blockSize)
            , _base(base)
            , _heap()
            , _id(0)
        {
            assert(!_base.empty());
            if (_base.back() != std::filesystem::path::preferred_separator)
                _base.push_back(std::filesystem::path::preferred_separator);
        }

        memory::memory_block MMapBlockAllocator::allocate_block()
        {
            uint64_t blockId = _id++;
            FilePointer blockPtr = std::make_unique<File>(
                _base + std::to_string(blockId) + ".mmap",
                _bsize);
            void *key = reinterpret_cast<void *>(blockPtr->data());

            typename Heap::accessor p;
            _heap.insert(p, key);
            p->second.swap(blockPtr);

            return { key, _bsize };
        }

        void MMapBlockAllocator::deallocate_block(memory::memory_block b)
        { _heap.erase(b.memory); }

        size_t MMapBlockAllocator::next_block_size() const noexcept
        { return _bsize; }

        DynamicAllocator makeDynamicAllocator(size_t maxRAMAlloc,
                                              size_t maxDiskAllocL2,
                                              std::string_view baseDir)
        {
            size_t maxAlloc,
                   blockSize,
                   maxRAMAllocL2 = ceil_log2(maxRAMAlloc);
            if (maxRAMAllocL2 < maxDiskAllocL2)
            {
                maxAlloc = (size_t)1u << maxDiskAllocL2;
                blockSize = maxDiskAllocL2 * maxAlloc;
            }
            else
            {
                maxAlloc = (size_t)1u << maxRAMAllocL2;
                blockSize = maxRAMAllocL2 * maxAlloc;
            }

            double mb = static_cast<double>(blockSize) / 1e6;
            LOG(INFO) << "created dynamic allocator with block size " << mb << "Mb";

            return DynamicAllocator {
                CappedAllocator(TbbMemoryPool(maxRAMAlloc, blockSize), maxRAMAlloc),
                memory::thread_safe_allocator<MMapMemoryPool, std::mutex> {
                    MMapMemoryPool(maxAlloc, blockSize, baseDir)
                }
            };
        }

    } // namespace alloc

    namespace detail
    {
        
        bool SharedGate::tryObtain()
        {
            Status expected = OPEN;
            return status.compare_exchange_strong(expected, CLOSED);
        }

        void SharedGate::release()
        {
            status = OPEN;
            status.notify_all();
        }

        void SharedGate::wait()
        {
            status.wait(CLOSED);
        }

    } // namespace detail
    
    RankClustering::RankClustering(size_t numColours)
        : coloursObserved(numColours)
        , colours(numColours + 1u, 1u)
        , ranks(numColours)
    {}

    size_t RankClustering::rank(size_t colour_) {
        if (coloursObserved.trySet(colour_)) {
            uint64_t rnk = colours.allocNextPosition();
            colours[rnk] = colour_;
            ranks[colour_].store(rnk);
            return rnk;
        } else {
            uint64_t rnk = 0;
            while (!rnk)
                rnk = ranks[colour_].load();
            return rnk;
        }
    }

    size_t RankClustering::colour(size_t rank_) const
    { return colours[rank_]; }

} // namespace pmap
