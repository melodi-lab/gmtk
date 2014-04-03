
/*
 * popen_err.c
 * 
 * A version of popen() that also makes stderr available
 * so that warnings/error messages can be prefixed.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include "popen_err.h"

#if HAVE_WORKING_FORK && HAVE_WAIT && HAVE_FDOPEN && HAVE_DUP2 && HAVE_EXECVP && HAVE_PIPE
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

/*
 * Normal popen() only gives us access to the forked process' stdout.
 * To address https://j.ee.washington.edu/trac/gmtk/ticket/389 we need
 * the process' stderr as well. So popen_err() gives us access to both.
 * Any lines of output on the process' stderr are prefixed by prefix.
 */
 
FILE *
popen_err(char const *command, char const *type, char const *prefix) {
  int pipefd[2], errfd[2];
  pid_t outid, errid;

  int argc;
  char const **argv = NULL;

  /* Create a pipe for the spawned process' stdout and another for stderr */

  if (pipe(pipefd) == -1) {
    perror("pipe");
    exit(EXIT_FAILURE);
  }
  if (pipe(errfd) == -1) {
    perror("pipe");
    exit(EXIT_FAILURE);
  }

  /* fork off a child process to run the command */

  outid = fork();
  if (-1 == outid) {
    perror("fork");
    exit(EXIT_FAILURE);
  }

  if (0 == outid) { /* child process - run command */

    /* close unneeded read end of pipes */
    if (close(pipefd[0]) || close(errfd[0])) {
      perror("closing pipe");
      _exit(EXIT_FAILURE);
    }
    /* make the pipes our stdout and stderr */
    if (dup2(pipefd[1], 1) == -1 || dup2(errfd[1],2) == -1) {
      perror("dup2");
      _exit(EXIT_FAILURE);
    }
    /* now run command */
#if 1 
    argc = countArgs(command);
    argv = (char const **)calloc(argc, sizeof(char const*));
    getArgs(command, argv);
    execvp(command, argv);
    perror("execvp");
    _exit(EXIT_FAILURE); /* only get here if execvp fails */ 

#else
    int i;
    for (i=0; i < 10; i+=1) {
      write(1, "12345\n", 6);
      write(2, "67890\n", 6);
    }

    /* all done - cleanup */
    if (close(pipefd[1]) || close(errfd[1])) {
      perror("closing pipe");
      _exit(EXIT_FAILURE);
    }
    _exit(EXIT_SUCCESS);
#endif
  } else { /* parent process - echo stderr & return read end of stdout pipe */

    /* close unneeded write end of pipes */
    if (close(pipefd[1]) || close(errfd[1])) {
      perror("closing pipe");
      exit(EXIT_FAILURE);
    }
    /* It would be painful for clients to have to read both stdout & stderr
     * (where do you read from them relative to each other?). So, we fork
     * off another process to asynchronously read the forked process' stderr
     * and echo it prefixed by prefix.
     */
    errid = fork();
    if (-1 == errid) {
      perror("fork");
      _exit(EXIT_FAILURE);
    }
    if (0 == errid) { /* child process - echo stderr with prefix */
      char buf[POPEN_BUF_SIZE+1];
      FILE *err = fdopen(errfd[0], "r");
      if (!err) {
        perror("fdopen");
        _exit(EXIT_FAILURE);
      }
      while (!feof(err) && !ferror(err)) {
        fgets(buf, POPEN_BUF_SIZE, err);
        fprintf(stderr,"%s%s", prefix, buf);
      }
      /* all done - cleanup */
      if (fclose(err)) {
        perror("fclose");
        _exit(EXIT_FAILURE);
      }
      _exit(EXIT_SUCCESS);
    }

    /* parent process - return the read end of the forked process' stdout
     * pipe to the client
     */ 
    return fdopen(pipefd[0], "r");
  }

  return NULL;
}

int 
pclose_err(FILE *stream) {
  int stat1, stat2;
  wait(&stat1); wait(&stat2); /* see how our children fared */
  /* The second child process (the stderr echoer) should have closed
   * the forked process' stderr fd. The only open file left should be
   * the forked process' stdout fd, so we close that.
   */
  return fclose(stream);
}

#else

FILE *
popen_err(char const *command, char const *type, char const *prefix) {
  return popen(command, type);
}

int 
pclose_err(FILE *stream) {
  return pclose(stream);
}

#endif
