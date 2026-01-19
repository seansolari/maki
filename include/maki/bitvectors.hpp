#pragma once

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <sdsl/vectors.hpp>
#include <maki/utils.hpp>

#include <mutex>

namespace fs = std::filesystem;

// https://en.cppreference.com/w/cpp/atomic/atomic_flag.html
class SpinLock
{
    std::atomic_flag _m{};
public:
    void lock() noexcept;
    bool try_lock() noexcept;
    void unlock() noexcept;
};

namespace bv
{
    
    namespace detail
    {
        
        struct DeadlockException : std::runtime_error {
            using std::runtime_error::runtime_error;
        };

        template <typename Mutex>
        class HeavyMutexManager {
        protected:
            std::vector<Mutex> mutexes;

        public:
            HeavyMutexManager() =default;
            HeavyMutexManager(size_t numLocks) : mutexes(numLocks) {}

            inline void realloc(size_t numLocks) { std::vector<Mutex>(numLocks).swap(mutexes); }
            inline size_t locks() const noexcept { return mutexes.size(); }
            inline size_t capacity() const noexcept { return mutexes.size(); }
            inline bool tryLock(size_t lockIndex) { return mutexes[lockIndex].try_lock(); }
            inline void lock(size_t lockIndex) { mutexes[lockIndex].lock(); }
            inline void unlock(size_t lockIndex) { mutexes[lockIndex].unlock(); }
        };

        using StdMutexManager = HeavyMutexManager<std::mutex>;
        using SpinMutexManager = HeavyMutexManager<SpinLock>;

        class MutexManager {
            size_t numWords,
                   numLocks;
            static constexpr size_t BITS_PER_WORD = 8;
            std::unique_ptr<std::atomic_uint8_t[]> _locks;

        public:
            MutexManager() =default;
            MutexManager(size_t numLocks);

            void realloc(size_t numLocks_);
            inline size_t locks() const noexcept { return numLocks; }
            inline size_t capacity() const noexcept { return 8 * numWords; }
            inline std::atomic_uint8_t* raw() const { return _locks.get(); }
            bool tryLock(size_t);
            void lock(size_t);
            void unlock(size_t);
            bool probe() const;
        };

    } // namespace detail
    
    namespace threaded
    {

        template <typename MutexManager>
        class LockedRegionManager {
        protected:
            size_t       regionSize;
            MutexManager mutexes;

        public:
            LockedRegionManager(size_t elems_, size_t regionSize_)
                : regionSize(regionSize_)
                , mutexes((elems_ + regionSize_ - 1) / regionSize_)
            {}

            LockedRegionManager() =default;

            void realloc(size_t elems_, size_t regionSize_) {
                regionSize = regionSize_;
                mutexes.realloc((elems_ + regionSize_ - 1) / regionSize_);
            }

        public:
            class Accessor {
                MutexManager &lockMgr;
                size_t       currentRegion,
                             regionSize;
                bool         hasLock;

            public:
                Accessor(MutexManager &mxs_, size_t regionSize)
                    : lockMgr(mxs_)
                    , currentRegion(-1)
                    , regionSize(regionSize)
                    , hasLock(false)
                {}

                Accessor(LockedRegionManager &mgr) : Accessor(mgr.mutexes, mgr.regionSize) {}

                static constexpr struct shallow_copy_tag_t{} shallow_copy_tag{};

                Accessor(Accessor const &rhs, shallow_copy_tag_t)
                    : lockMgr(rhs.lockMgr)
                    , currentRegion(-1)
                    , regionSize(rhs.regionSize)
                    , hasLock(false)
                {}

                Accessor(Accessor const&) =delete;
                Accessor& operator=(Accessor const&) =delete;

                ~Accessor() {
                    releaseLock();
                }

                void access(size_t index) {
                    size_t region = index / regionSize;
                    if (region != currentRegion) {
                        releaseLock();
                        lockMgr.lock(region);
                        currentRegion = region;
                        hasLock = true;
                    }
                }

                void releaseLock() {
                    if (hasLock) {
                        lockMgr.unlock(currentRegion);
                        currentRegion = -1;
                        hasLock = false;
                    }
                }
            };
        };

        using SpinRegionManager = LockedRegionManager<detail::SpinMutexManager>;

        template <uint8_t width, typename MutexManager>
        class ManagedIntVector : public LockedRegionManager<MutexManager> {
        private:
            static_assert(width > 0);

        protected:
            using ParentAccessor = typename LockedRegionManager<MutexManager>::Accessor;
            using BaseVector = sdsl::int_vector<width>;
            BaseVector _vec;

        public:
            using Rank1Type   = typename BaseVector::rank_1_type;
            using Rank0Type   = typename BaseVector::rank_0_type;
            using Select1Type = typename BaseVector::select_1_type;
            using Select0Type = typename BaseVector::select_0_type;

            inline void setToValue(uint64_t v)
            { sdsl::util::set_to_value(_vec, v); }

            template <typename T>
            inline void initSupport(T &obj)
            { sdsl::util::init_support(obj, &_vec); }

            inline auto sizeInBytes() const
            { return sdsl::size_in_bytes(_vec); }

        public:
            class ConstAccessor;
            class Accessor;

            ManagedIntVector() =default;

            ManagedIntVector(size_t size, size_t grainsize)
                : LockedRegionManager<MutexManager>::LockedRegionManager(size, grainsize)
                , _vec(size, 0)
            {
                assert((grainsize * width) % (8 * sizeof(typename BaseVector::value_type)) == 0);
            }

            ManagedIntVector(BaseVector const &data, size_t grainsize)
                : LockedRegionManager<MutexManager>::LockedRegionManager(data.size(), grainsize)
                , _vec(data)
            {
                assert((grainsize * width) % (8 * sizeof(typename BaseVector::value_type)) == 0);
            }

            ManagedIntVector(BaseVector &&data, size_t grainsize)
                : LockedRegionManager<MutexManager>::LockedRegionManager(data.size(), grainsize)
                , _vec(std::move(data))
            {
                assert((grainsize * width) % (8 * sizeof(typename BaseVector::value_type)) == 0);
            }

            ManagedIntVector(ManagedIntVector const &rhs)
                : LockedRegionManager<MutexManager>::LockedRegionManager(rhs.size(), rhs.regionSize)
                , _vec(rhs._vec)
            {}

            inline void saveToDisk(fs::path const &file) const
            { fileutils::store_to_file(_vec, file); }

            void loadFromDisk(fs::path const &file, size_t grainsize) {
                fileutils::load_from_file(_vec, file);
                this->realloc(_vec.size(), grainsize);
            }

        public:
            class ConstAccessor {
                BaseVector const &vec;

            public:
                ConstAccessor(BaseVector const &p) : vec(p) {}
                ConstAccessor(ManagedIntVector const &bv) : ConstAccessor(bv._vec) {}
                inline operator BaseVector const&() const noexcept { return vec; }
                inline auto operator[](size_t const &i) const { return vec[i]; }
                inline auto begin() const { return vec.begin(); }
                inline auto end() const { return vec.end(); }
                inline auto cbegin() const { return vec.cbegin(); }
                inline auto cend() const { return vec.cend(); }
            };

            class Accessor : public ParentAccessor {
                BaseVector &vec;
            
            public:
                Accessor(BaseVector &v, MutexManager &mxs_, size_t regionSize)
                    : ParentAccessor(mxs_, regionSize)
                    , vec(v)
                {}

                Accessor(ManagedIntVector &vec)
                    : Accessor(vec._vec, vec.mutexes, vec.regionSize)
                {}

                Accessor(Accessor const &rhs, typename ParentAccessor::shallow_copy_tag_t)
                    : ParentAccessor(rhs, ParentAccessor::shallow_copy_tag)
                    , vec(rhs.vec)
                {}

                Accessor(Accessor const&) =delete;
                Accessor& operator=(Accessor const&) =delete;

                inline auto operator[](size_t index) {
                    this->access(index);
                    return vec[index];
                }
            };

        public:
            inline void operator|=(ManagedIntVector const &rhs)
            { _vec |= rhs._vec; }

            inline BaseVector&& data() &&
            { return std::move(_vec); }

        public:
            inline auto size() const noexcept
            { return _vec.size(); }

            inline void resize(size_t size) {
                _vec.resize(size);
                this->realloc(size, this->regionSize);
            }
        };

        template <uint8_t width>
        using IntVector = ManagedIntVector<width,detail::SpinMutexManager>;

        class ConcurrentBitVector
        {
            size_t _size;
            std::vector<copyable_atomic_uint64_t> _data;

        public:
            ConcurrentBitVector(size_t nbits_);

            /** trySet
             * Try to set bit `bit`. If it is not set, then set the bit. Returns whether the bit was set.
             */
            bool trySet(size_t bit);

            inline std::vector<copyable_atomic_uint64_t> const &raw() const noexcept
            { return _data; }
        };

    } // namespace threaded
    
} // namespace bv
