
#include "maki/core/utils/locks.hpp"

MutexManager::MutexManager(size_t numLocks)
    : numWords((numLocks + BITS_PER_WORD - 1) / BITS_PER_WORD),
      numLocks(numLocks),
      _locks(std::make_unique<std::atomic_uint8_t[]>(numWords)) {}

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
    if (_locks[wordIndex].compare_exchange_weak(expected, expected | mask,
                                                std::memory_order_acquire,
                                                std::memory_order_relaxed)) {
      return true;
    }
  }
  return false;
}

void MutexManager::lock(size_t lockIndex) {
  while (!tryLock(lockIndex)) {
  }
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
    if ((p + i)->load() != 0)
      return true;
  return false;
}

void LockedRegionManager::realloc(size_t elems_, size_t regionSize_) {
  regionSize = regionSize_;
  mutexes.realloc((elems_ + regionSize_ - 1) / regionSize_);
}

void LockedRegionManager::Accessor::access(size_t index) {
  size_t region = index / regionSize;
  if (region != currentRegion) {
    releaseLock();
    lockMgr.lock(region);
    currentRegion = region;
    hasLock = true;
  }
}

void LockedRegionManager::Accessor::releaseLock() {
  if (hasLock) {
    lockMgr.unlock(currentRegion);
    currentRegion = -1;
    hasLock = false;
  }
}
