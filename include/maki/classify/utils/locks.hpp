
#pragma once
#include <atomic>
#include <memory>

class MutexManager {
  std::size_t numWords, numLocks;
  static constexpr std::size_t BITS_PER_WORD = 8;
  std::unique_ptr<std::atomic_uint8_t[]> _locks;

public:
  MutexManager() = default;
  MutexManager(std::size_t numLocks);

  void realloc(std::size_t numLocks_);
  inline std::size_t locks() const noexcept { return numLocks; }
  inline std::size_t capacity() const noexcept { return 8 * numWords; }
  inline std::atomic_uint8_t *raw() const { return _locks.get(); }
  bool tryLock(std::size_t);
  void lock(std::size_t);
  void unlock(std::size_t);
  bool probe() const;
};

class LockedRegionManager {
protected:
  size_t regionSize;
  MutexManager mutexes;

public:
  LockedRegionManager(size_t elems_, size_t regionSize_)
      : regionSize(regionSize_),
        mutexes((elems_ + regionSize_ - 1) / regionSize_) {}
  LockedRegionManager() = default;
  
  void realloc(size_t elems_, size_t regionSize_);

public:
  class Accessor {
  public:
    Accessor(MutexManager &mxs_, size_t regionSize)
        : lockMgr(mxs_), currentRegion(-1), regionSize(regionSize),
          hasLock(false) {}
    Accessor(LockedRegionManager &mgr)
        : Accessor(mgr.mutexes, mgr.regionSize) {}
    static constexpr struct shallow_copy_tag_t {
    } shallow_copy_tag{};
    Accessor(Accessor const &rhs, shallow_copy_tag_t)
        : lockMgr(rhs.lockMgr), currentRegion(-1), regionSize(rhs.regionSize),
          hasLock(false) {}
    Accessor(Accessor const &) = delete;
    Accessor &operator=(Accessor const &) = delete;
    ~Accessor() { releaseLock(); }

    void access(size_t index);
    void releaseLock();

  protected:
    MutexManager &lockMgr;
    size_t currentRegion, regionSize;
    bool hasLock;
  };
};