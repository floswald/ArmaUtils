
#ifndef _ArmaUtils_HANDLER_H
#define _ArmaUtils_HANDLER_H

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>



inline void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

#endif

