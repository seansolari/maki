#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>

/* ===================== CONFIG ===================== */

// Compile-time minimum log level
#ifndef LOG_COMPILED_LEVEL
#define LOG_COMPILED_LEVEL 1 // 0=Debug, 1=Info, 2=Warn, 3=Error
#endif

// Async queue size limit
#ifndef LOG_QUEUE_LIMIT
#define LOG_QUEUE_LIMIT 8192
#endif

// Use UTC timestamps
#ifndef LOG_USE_UTC
#define LOG_USE_UTC 1
#endif

// Log rotation (bytes, 0 = disabled)
#ifndef LOG_ROTATE_BYTES
#define LOG_ROTATE_BYTES (5 * 1024 * 1024)
#endif

/* ================================================== */

class Logger {
public:
  enum class Level { Debug = 0, Info, Warning, Error };

  static Logger &instance() {
    static Logger inst;
    return inst;
  }

  Logger(const Logger &) = delete;
  Logger &operator=(const Logger &) = delete;

  void enable_file(const std::string &path) {
    std::lock_guard<std::mutex> lock(file_mutex_);
    file_path_ = path;
    file_.open(path, std::ios::app);
  }

  bool enabled(Level level) const {
    return static_cast<int>(level) >= LOG_COMPILED_LEVEL;
  }

  void enqueue(Level level, std::string message) {
    if (!enabled(level))
      return;

    std::unique_lock<std::mutex> lock(queue_mutex_);
    queue_cv_.wait(lock, [&] { return queue_.size() < LOG_QUEUE_LIMIT; });

    queue_.push({level, std::move(message)});
    lock.unlock();
    queue_cv_.notify_one();
  }

  class Stream {
  public:
    Stream(Level level) : level_(level) {}
    ~Stream() { Logger::instance().enqueue(level_, stream_.str()); }

    template <typename T> Stream &operator<<(const T &value) {
      stream_ << value;
      return *this;
    }

  private:
    Level level_;
    std::ostringstream stream_;
  };

private:
  struct Entry {
    Level level;
    std::string message;
  };

  Logger() : stop_(false) {
    worker_ = std::thread([this] { worker_loop(); });
  }

  ~Logger() {
    stop_ = true;
    queue_cv_.notify_all();
    if (worker_.joinable())
      worker_.join();
  }

  void worker_loop() {
    while (!stop_) {
      std::unique_lock<std::mutex> lock(queue_mutex_);
      queue_cv_.wait(lock, [&] { return stop_ || !queue_.empty(); });

      while (!queue_.empty()) {
        Entry e = std::move(queue_.front());
        queue_.pop();
        lock.unlock();

        write(e);

        lock.lock();
        queue_cv_.notify_one();
      }
    }
  }

  void write(const Entry &e) {
    const std::string line =
        timestamp() + " [" + level_str(e.level) + "] " + e.message + "\n";

    std::ostream &console =
        (e.level == Level::Error || e.level == Level::Warning) ? std::cerr
                                                               : std::cout;

    console << line;
    console.flush();

    std::lock_guard<std::mutex> lock(file_mutex_);
    if (file_.is_open()) {
      rotate_if_needed();
      file_ << line;
      file_.flush();
    }
  }

  std::string timestamp() const {
    using clock = std::chrono::system_clock;
    auto now = clock::now();
    std::time_t t = clock::to_time_t(now);

    std::tm tm{};
#if LOG_USE_UTC
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
#else
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
  }

  std::string level_str(Level lvl) const {
    switch (lvl) {
    case Level::Debug:
      return "DEBUG";
    case Level::Info:
      return "INFO";
    case Level::Warning:
      return "WARN";
    case Level::Error:
      return "ERROR";
    }
    return "UNKNOWN";
  }

  void rotate_if_needed() {
#if LOG_ROTATE_BYTES > 0
    if (file_.tellp() >= LOG_ROTATE_BYTES) {
      file_.close();
      std::string rotated = file_path_ + ".old";
      std::remove(rotated.c_str());
      std::rename(file_path_.c_str(), rotated.c_str());
      file_.open(file_path_, std::ios::trunc);
    }
#endif
  }

private:
  std::queue<Entry> queue_;
  std::mutex queue_mutex_;
  std::condition_variable queue_cv_;
  std::atomic<bool> stop_;
  std::thread worker_;

  std::ofstream file_;
  std::string file_path_;
  std::mutex file_mutex_;
};

/* ===================== MACROS ===================== */

#define LOG_DEBUG()                                                            \
  if (!Logger::instance().enabled(Logger::Level::Debug))                       \
    ;                                                                          \
  else                                                                         \
    Logger::Stream(Logger::Level::Debug)

#define LOG_INFO()                                                             \
  if (!Logger::instance().enabled(Logger::Level::Info))                        \
    ;                                                                          \
  else                                                                         \
    Logger::Stream(Logger::Level::Info)

#define LOG_WARN()                                                             \
  if (!Logger::instance().enabled(Logger::Level::Warning))                     \
    ;                                                                          \
  else                                                                         \
    Logger::Stream(Logger::Level::Warning)

#define LOG_ERROR()                                                            \
  if (!Logger::instance().enabled(Logger::Level::Error))                       \
    ;                                                                          \
  else                                                                         \
    Logger::Stream(Logger::Level::Error)

/* ================================================== */
