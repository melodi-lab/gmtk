/*
 * "Copyright 2001, International Business Machines Corporation and University
 * of Washington. All Rights Reserved
 *
 *    Written by Geoffrey Zweig and Jeff Bilmes
 *
 * NO WARRANTY
 * THE PROGRAM IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED INCLUDING, WITHOUT
 * LIMITATION, ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT,
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Each Recipient is
 * solely responsible for determining the appropriateness of using the Program
 * and assumes all risks associated with such use, including but not limited
 * to the risks and costs of program errors, compliance with applicable laws,
 * damage to or loss of data, programs or equipment, and unavailability or
 * interruption of operations.

 * DISCLAIMER OF LIABILITY
 * THE UNIVERSITY OF WASHINGTON, INTERNATIONAL BUSINESS MACHINES CORPORATION,
 * GEOFFREY ZWEIG AND JEFF BILMES SHALL NOT HAVE ANY LIABILITY FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING WITHOUT LIMITATION LOST PROFITS), HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  OF
 * THE PROGRAM, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."
*/

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <vector>
#include <string>
#include <map>

/* A class to read argument lists. 
   The main structure is like this:
   If the argument "-argsFile" is specified, 
       then parameters are read from the file
   If the argument "-argsFile" is specified more than once, 
       then parameters are read from each file
   (-argsFile must be the name of the file-reading flag)
   The remaining arguments are read and parsed.
   (This allows file-parameters to be "overrided" by the command line.)
   If an argument is specified multiple times, then the last
   specification is used, and a warning is printed.
   Unambiguous flag prefixes may be used instead of the full name.
   Checking:
       If a required argument is not specifed an error results.
       Some simple checking is done to verify that strings look like
       the types they are supplosed to represent. for example, if "count"
       is an integer, then "-count 1.0" will be identified as an error.
       The use of ambiguous prefixes causes an error.
       Recursive file specification signals an error.
*/        

/*
   Known bug: quoted strings with spaces will not be properly interpreted if
              read from a file.
*/
    
struct Argument
{
    string name;
    string type;
    bool required;
    void *storage;
    string remark;
    Argument(string n, string t, bool req, void *s, string r)
        {name=n; type=t; required=req; storage=s; remark=r;}
};

struct Argument_List
{
    vector<Argument> arguments;
    map<string, int> countFor;
    
    void add(string _name, bool _required, int *s, string desc);
    void add(string _name, bool _required, unsigned *s, string desc);
    void add(string _name, bool _required, float *s, string desc);
    void add(string _name, bool _required, double *s, string desc);
    void add(string _name, bool _required, bool *s, string desc);
    void add(string _name, bool _required, char **s, string desc);
 
    void parse(int count, char *args[]); 

    string fullName(string n);  
    // expands a name-prefix n into the full name
    // e.g. baseCaseT => baseCaseThreshold

    Argument *argumentNamed(string s);
    void getVal (Argument *arg, char *val);
    void dumpArgs();
    void check();
    void getArgVals(int count, char *args[]);
    void checkArgStructure(int count, char *args[]);
    void readFile(char *fn, int &c, char ** &av);
    void checkRedundancy(string n);
};

#endif
