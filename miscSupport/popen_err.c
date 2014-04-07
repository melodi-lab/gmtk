
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
#include <string.h>
#include <assert.h>

#include "popen_err.h"

#if HAVE_WORKING_FORK && HAVE_WAIT && HAVE_FDOPEN && HAVE_DUP2 &&	\
    HAVE_EXECVP && HAVE_PIPE && HAVE_STRCPY && HAVE_STRTOK_R && HAVE_STRDUP

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

static int
allSpace(char const *s) {
  for (; *s; s+=1) 
    if (*s != ' ') 
      return 0;
  return -1;
}

static void
parseArgs(char const *command, char *argv[]) {
  char *s, *arg, *last;
  int argc = 0;
  assert(command);
  assert(argv);
  s = strdup(command); /* no need to free, process will end */
  assert(s);
  /* turn line continuations into white space */
  for (last = strstr(s, "\\\n"); last; last = strstr(last, "\\\n")) {
    last[0]=' ';
    last[1]=' ';
  }
  /* skip white space */
  for ( ; *s == ' '; s+=1)
    ;
  arg = strtok_r(s, " ", &last);
  argv[argc++] = arg;
  do {
    if (argc >= POPEN_MAX_ARGC) {
      fprintf(stderr, "too many arguments to popen_err('%s'), there can be at most %d\n", 
	      command, POPEN_MAX_ARGC);
      _exit(EXIT_FAILURE);
    }
    arg = strtok_r(NULL, " ", &last);
    if (arg && allSpace(arg)) 
      continue;
    argv[argc++] = arg;
  } while (arg);
}


/*
 * Normal popen() only gives us access to the forked process' stdout.
 * To address https://j.ee.washington.edu/trac/gmtk/ticket/389 we need
 * the process' stderr as well. So popen_err() gives us access to both.
 * Any lines of output on the process' stderr are prefixed by prefix.
 */
 
FILE *
popen_err(char const *command, char const *type, char const *prefix) {
  int outfd[2], errfd[2];
  pid_t outid, errid;
  char *argv[POPEN_MAX_ARGC+1];

  /* Create a pipe for the spawned process' stdout and another for stderr */

  if (pipe(outfd) == -1) {
    perror("popen_err(): failed to create stdout pipe");
    return NULL;
  }
  if (pipe(errfd) == -1) {
    perror("popen_err(): failed to create stderr pipe");
    return NULL;
  }

  /* fork off a child process to run the command */

  outid = fork();
  if (-1 == outid) {
    perror("popen_err(): failed to create process to execute command");
    return NULL;
  }

  if (0 == outid) { /* child process - run command */

    /* close unneeded read end of pipes */
    if (close(outfd[0]) || close(errfd[0])) {
      perror("popen_err(): failed to close read end of pipe in command process");
      _exit(EXIT_FAILURE);
    }

    /* make the pipes our stdout and stderr */
    if (dup2(outfd[1], STDOUT_FILENO) == -1 || dup2(errfd[1], STDERR_FILENO) == -1) {
      perror("popen_err(): failed to connect command process stdout or stderr to pipe");
      _exit(EXIT_FAILURE);
    }

    /* now run command */
    parseArgs(command, argv); 
    execvp(argv[0], argv);
    /* only get here if execvp fails */
    perror("popen_err(): failed to execute command");
    _exit(EXIT_FAILURE);

  } else { /* parent process - echo stderr & return read end of stdout pipe */
    /* close unneeded write end of pipes */
    if (close(outfd[1]) || close(errfd[1])) {
      perror("popen_err(): failed to close write end of stdout or stderr pipe(s)");
      return NULL;
    }
    /* It would be painful for clients to have to read both stdout & stderr
     * (where do you read from them relative to each other?). So, we fork
     * off another process to asynchronously read the forked process' stderr
     * and echo it prefixed by prefix.
     */
    errid = fork();
    if (-1 == errid) {
      perror("popen_err(): failed to create stderr echo process");
      return NULL;
    }
    if (0 == errid) { /* child process - echo stderr with prefix */
      char buf[POPEN_BUF_SIZE+1];
      FILE *err = fdopen(errfd[0], "r");
      if (!err) {
        perror("popen_err(): openning stderr stream in echo process failed");
        _exit(EXIT_FAILURE);
      }
      while (1) {
        fgets(buf, POPEN_BUF_SIZE, err);
        if (feof(err) || ferror(err)) break;
        fprintf(stderr,"%s%s", prefix, buf);
      }
      /* all done - cleanup */
      if (fclose(err)) {
        perror("popen_err(): closing stderr stream in echo process failed");
        _exit(EXIT_FAILURE);
      }
      _exit(EXIT_SUCCESS);
    }
    /* parent process - return the read end of the forked process' stdout
     * pipe to the client
     */ 
    return fdopen(outfd[0], "r");
  }

  return NULL;
}

int 
pclose_err(FILE *stream) {
  int stat1, stat2, result;
  wait(&stat1); wait(&stat2); /* see how our children fared */
  if (stat1 || stat2) {
    fprintf(stderr,"pclose_err(): child process failed: %s\n", strerror(errno));
    return -1;
  }
  /*
   * The second child process (the stderr echoer) should have closed
   * the forked process' stderr fd. The only open file left should be
   * the forked process' stdout fd, so we close that.
   */
  result = fclose(stream);
  if (result) {
    perror("pclose_err(): closing stdout stream failed");
    return -1;
  }
  return 0;
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
