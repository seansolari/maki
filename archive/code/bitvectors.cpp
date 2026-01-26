#include <maki/bitvectors.hpp>

#include <glog/logging.h>

void SpinLock::lock() noexcept
{
    while (_m.test_and_set(std::memory_order_acquire))
        _m.wait(true, std::memory_order_relaxed);
}

bool SpinLock::try_lock() noexcept
{
    return !_m.test_and_set(std::memory_order_acquire);
}

void SpinLock::unlock() noexcept
{
    if (!_m.test(std::memory_order_acquire))
        LOG(INFO) << "ERROR: releasing flag when not acquired!";
    _m.clear(std::memory_order_release);
    _m.notify_one();
}

namespace bv
{

    namespace detail
    {
    
        MutexManager::MutexManager(size_t numLocks)
            : numWords((numLocks + BITS_PER_WORD - 1) / BITS_PER_WORD)
            , numLocks(numLocks)
            , _locks(std::make_unique<std::atomic_uint8_t[]>(numWords))
        {
        }

        void MutexManager::realloc(size_t numLocks_) {
            numWords = (numLocks_ + BITS_PER_WORD - 1) / BITS_PER_WORD;
            numLocks = numLocks_;
            std::make_unique<std::atomic_uint8_t[]>(numWords).swap(_locks);
        }

        bool MutexManager::tryLock(size_t lockIndex) {
            size_t wordIndex = lockIndex / BITS_PER_WORD,
                   bitOffset = lockIndex % BITS_PER_WORD;
            uint8_t mask = 1U << bitOffset,
                    expected = _locks[wordIndex].load(std::memory_order_relaxed);
            while (!(expected & mask)) {
                if (_locks[wordIndex].compare_exchange_weak(
                        expected, expected | mask,
                        std::memory_order_acquire,
                        std::memory_order_relaxed)) {
                    return true;
                }
            }
            return false;
        }

        void MutexManager::lock(size_t lockIndex) {
            while (!tryLock(lockIndex)) {}
        }

        void MutexManager::unlock(size_t lockIndex) {
            size_t wordIndex = lockIndex / BITS_PER_WORD,
                   bitOffset = lockIndex % BITS_PER_WORD;
            uint8_t mask = ~(1U << bitOffset);
            _locks[wordIndex].fetch_and(mask, std::memory_order_release);
        }

        bool MutexManager::probe() const {
            auto const *p = raw();
            for (size_t i = 0; i < numWords; ++i)
                if ((p+i)->load() != 0) return true;
            return false;
        }

    } // namespace detail

    namespace threaded
    {
        
        ConcurrentBitVector::ConcurrentBitVector(size_t nbits_)
            : _size((nbits_ + 63) / 64),
            _data(_size, 0u)
        {
        }

        bool ConcurrentBitVector::trySet(size_t bit) {
            size_t celli = bit / 64,
                   biti = bit % 64;

            assert(celli < _size);
            auto atom = _data[celli].ref();

            uint64_t current = atom.load(),
                     mask = 1ull << biti;
            bool swapped = false;
            while (((current & mask) == 0) && (!swapped))
                swapped = atom.compare_exchange_weak(current, current | mask);

            return swapped;
        }


    } // namespace threaded
    
} // namespace bv
