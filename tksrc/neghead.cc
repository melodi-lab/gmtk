
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

/*
 * This is a test program that provides the funcationality
 * of GNU head's negative arguments: "head -n -3 file" prints
 * all but the last 3 lines of file. Unfortunately, non-GNU
 * head implementations don't support this feature (in 
 * particular the default OSX head does not.
 *
 */

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "neghead n" << endl;
    return 1;
  }
  int n = atoi(argv[1]);
  if (n < 1) {
    cerr << "first argument must be > 0" << endl;
    return 1;
  }
  string *lines = new string[n];
  if (lines == NULL) {
    cerr << "unable to allocate " << n << " strings" << endl;
    return 1;
  }

  // read first n lines
  int j;
  for (j=0; j < n && cin.good(); j+=1)
    getline(cin, lines[j]);
  // now print the earliest line read and replace it 
  // with the next line read
  for (j=0; cin.good(); j = (j + 1) % n) {
    string line;
    getline(cin,line);
    if (cin.good()) {
      cout << lines[j] << endl;
      lines[j] = line;
    }
  }
  
  return 0;
}
