#pragma once

// Small POSIX process/file facade for SingularInput.cpp.
// These POSIX APIs have no standard C++ equivalent, so the headers are kept
// behind this single facade header.

#include <fcntl.h>     // open, O_RDONLY
#include <sys/types.h> // pid_t
#include <sys/wait.h>  // waitpid, WIFEXITED, WEXITSTATUS
#include <unistd.h>    // close, dup2, execvp, fork, mkstemp, STDIN_FILENO

namespace palp {

/// File descriptor wrapper for small POSIX file descriptors (close on scope
/// exit).
class PosixFd {
  int fd_;

public:
  explicit PosixFd(int fd) : fd_(fd) {}
  PosixFd(const PosixFd &) = delete;
  PosixFd &operator=(const PosixFd &) = delete;
  PosixFd(PosixFd &&other) noexcept : fd_(other.fd_) { other.fd_ = -1; }
  PosixFd &operator=(PosixFd &&other) noexcept {
    if (this != &other) {
      close();
      fd_ = other.fd_;
      other.fd_ = -1;
    }
    return *this;
  }
  ~PosixFd() { close(); }

  [[nodiscard]] int get() const noexcept { return fd_; }
  [[nodiscard]] bool valid() const noexcept { return fd_ != -1; }

  int release() noexcept {
    int tmp = fd_;
    fd_ = -1;
    return tmp;
  }

  void close() noexcept {
    if (fd_ != -1) {
      ::close(fd_);
      fd_ = -1;
    }
  }
};

/// Spawn a process replacing stdin with `scriptPath`, wait, and return exit
/// status. Returns true on clean exit(0); returns false on fork/exec/wait error
/// and prints `perror` for the caller.
[[nodiscard]] inline bool spawnSingular(const char *scriptPath,
                                        const char *const *argv) {
  pid_t pid = fork();
  if (pid == -1) {
    return false;
  }
  if (pid == 0) {
    // Child
    int fd = ::open(scriptPath, O_RDONLY);
    if (fd == -1 || ::dup2(fd, STDIN_FILENO) == -1 || ::close(fd) == -1) {
      _exit(1);
    }
    ::execvp(argv[0], const_cast<char *const *>(argv));
    _exit(127); // execvp failed
  }

  int status;
  if (waitpid(pid, &status, 0) == -1) {
    return false;
  }
  return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

} // namespace palp
