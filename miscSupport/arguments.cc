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

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <stdlib.h>
#include "arguments.h"

void Argument_List::checkRedundancy(string n)
{
    // make sure the argument isn't listed twice
    for (unsigned i=0; i<arguments.size(); i++)
        if (arguments[i].name == n)
        {
            cout << "Error: " << n << " argument added redundantly\n";
            exit(1);
        }
}

void Argument_List::add(string _name, bool _required, int *s, string desc)
{
    checkRedundancy(_name = "-" + _name);
    arguments.push_back(Argument(_name,"int",_required,s,desc));
}

void Argument_List::add(string _name, bool _required, unsigned *s, string desc)
{
    checkRedundancy(_name = "-" + _name);
    arguments.push_back(Argument(_name,"unsigned",_required,s,desc));
}

void Argument_List::add(string _name, bool _required, float *s, string desc)
{
    checkRedundancy(_name = "-" + _name);
    arguments.push_back(Argument(_name,"float",_required,s,desc));
}

void Argument_List::add(string _name, bool _required, double *s, string desc)
{
    checkRedundancy(_name = "-" + _name);
    arguments.push_back(Argument(_name,"double",_required,s,desc));
}

void Argument_List::add(string _name, bool _required, char **s, string desc)
{
    checkRedundancy(_name = "-" + _name);
    arguments.push_back(Argument(_name,"string",_required,s,desc));
}

void Argument_List::add(string _name, bool _required, bool *s, string desc)
{
    checkRedundancy(_name = "-" + _name);
    arguments.push_back(Argument(_name,"bool",_required,s,desc));
}

void Argument_List::check()
{
    // check that all required arguments have been read
    for (unsigned i=0; i<arguments.size(); i++)
        if (arguments[i].required == true)
            if (countFor[arguments[i].name] == 0)
            {
                cout << "Error: required argument named " << arguments[i].name
                     << " not specified\n";
                exit(1);
            }

    for (unsigned i=0; i<arguments.size(); i++)
        if (countFor[arguments[i].name]>1)
            cout << "Warning: argument " << arguments[i].name 
                 << " was specified " << countFor[arguments[i].name]
                 << " times; using the last\n";
}

Argument *Argument_List::argumentNamed(string s)
{
    for (unsigned i=0; i<arguments.size(); i++)
        if (arguments[i].name == s)
            return &arguments[i];
}

string Argument_List::fullName(string n)
{
    string fn = "";
    for (unsigned i=0; i<arguments.size(); i++)
        if (arguments[i].name.find(n) == 0) // n is a prefix
        {
            if (fn != "")
            {
                cout << "Error: " << n << " is a prefix of more than one argument\n";
                cout << "Legal arguments are: " << endl;
                dumpArgs();
                exit(1);
            }
            fn = arguments[i].name;
        }
    if (fn=="")
    {
        cout << "Error: unknown argument " << n << endl;
        cout << "Legal arguments are: " << endl;
        dumpArgs();
        exit(1);
    }
    return fn;
}

void verifyInt(char *s)
{
    for (char *p=s; *p!='\0'; p++)
        if (!isdigit(*p) && *p!='-')
        {
            cout << "Error: " << s << " does not look like an integer\n";
            exit(1);
        }
}

void verifyFloat(char *s)
{
    for (char *p=s; *p!='\0'; p++)
        if (!isdigit(*p) && *p!='-' && *p!='e' && *p!='.')
        {
            cout << "Error: " << s << " does not look like a float\n";
            exit(1);
        }
}

void Argument_List::getVal(Argument *arg, char *val)
{
    // interpret the string according to the argument type
    // more types can be added by listing them
    if (arg->type == "int")
    {
        verifyInt(val); 
        int v = atoi(val);
        *((int *)arg->storage) = v;
    }
    else if (arg->type=="unsigned")
    {
        verifyInt(val); 
        unsigned v = atoi(val);
        *((unsigned *)arg->storage) = v;
    }
    else if (arg->type=="float")
    {
        verifyFloat(val);
        float v = atof(val);
        *((float *)arg->storage) = v;
    }
    else if (arg->type=="double")
    {
        verifyFloat(val);
        double v = atof(val);
        *((double *)arg->storage) = v;
    }
    else if (arg->type=="string")
    {
        *((char **)arg->storage) = val;
    }
    else if (arg->type=="bool")
    {
        bool v;
        string vl(val);
        if (vl=="true" || vl=="True" || vl=="t" || vl=="T" || vl=="1")
            v = true;
        else if (vl=="false" || vl=="False" || vl=="f" || vl=="F" || vl=="0")
            v = false;
        else
        {
            cout << "Error: unknown value for boolean argument: " << val<<endl;
            cout << "Must be true;True;t;T;1;false;False;f;F;0\n";
            exit(1);
        }
        *((bool *)arg->storage) = v;
    }
    countFor[arg->name]++;
}

void Argument_List::dumpArgs()
{
    for (unsigned i=0; i<arguments.size(); i++)
    {
        cout << arguments[i].name << " " 
             << "(" << arguments[i].type;
        if (arguments[i].required)
            cout << ", required): ";
        else
        {
            cout << " [";
            if (arguments[i].type == "int")
                cout << *((int *)arguments[i].storage);
            else if (arguments[i].type=="unsigned")
                cout << *((unsigned *)arguments[i].storage);
            else if (arguments[i].type=="float")
                cout << *((float *)arguments[i].storage);
            else if (arguments[i].type=="double")
                cout << *((double *)arguments[i].storage);
            else if (arguments[i].type=="string")
                cout << *((char **)arguments[i].storage);
            else if (arguments[i].type=="bool")
                cout << *((bool *)arguments[i].storage);
            cout << "]): ";
        }
        cout << endl << "\t" << arguments[i].remark << endl;
    }
}
 
void Argument_List::checkArgStructure(int count, char *args[])
{
    if ((count % 2) == 0) 
    {
        cout << "Error: Argument pattern must be: prog_name (-arg_name arg_value)+\n";
        exit(1);
    }

    for (int i=1; i<count; i+=2)
        if (*args[i] != '-')
        {
            cout << "Error: argument names must begin with a -\n";
            cout << "See " << args[i] << endl;
            exit(1);
        }
}

void Argument_List::getArgVals(int count, char *args[])
{
    for (int i=1; i<count; i+=2)    
    {
        // get the name of the argument 
        string s = fullName(args[i]);
        
        if (s == "-argsFile")
            continue;

        // get the argument structure
        Argument *arg = argumentNamed(s);

        // get the argument value and shovel it into the right memory location
        getVal(arg, args[i+1]); 
    }
}

void Argument_List::readFile(char *fn, int &c, char ** &av)
{
    ifstream in(fn);
    if (!in) {cout << "Error: Unable to open " << fn << endl; exit(1);}
    cout << "Reading parameters from " << fn << endl;

    vector<string> arg_name;
    vector<string> arg_val;
    while (!in.eof())
    {
        string an, v;
        in >> an >> v >> ws;
        arg_name.push_back(an);
        arg_val.push_back(v);
    }

    if (arg_name.size() != arg_val.size())
    {
        cout << "Error: mismatched number of argument names and values in "
             << fn << endl;
        exit(1);
    }

    for (unsigned i=0; i<arg_name.size(); i++)
        if (fullName(arg_name[i]) == "-argsFile")
        {
            cout << "Error: Recursive file inclusion is disallowed. See "
                 << fn << endl;
            exit(1);
        } 

    // start numbering from 1 to match main(argc, argv[])
    c = 2*arg_name.size()+1;
    av = new char *[c];
    // move to the char *[]  format for consistency with args from main program
    for (unsigned i=0; i<arg_val.size(); i++)
    {
        int p = 2*i;
        av[p+1] = new char[arg_name[i].size()+1];
        strcpy(av[p+1], arg_name[i].c_str());
        av[p+2] = new char[arg_val[i].length()+1];
        strcpy(av[p+2], arg_val[i].c_str());
    }
   in.close();
}

void Argument_List::parse(int count, char *args[])
{
    if (count==1)
    {
        cout << "Argument usage: " << endl;
        dumpArgs();
        exit(1);
    }

    checkArgStructure(count, args);

    // first read argument values from any files
    for (int i=1; i<count; i+=2)
    {
        string s = fullName(args[i]);  
        if (s == "-argsFile")  // parse the arguments in the file
        {
            char **av;
            int ac;
            readFile(args[i+1], ac, av);
            checkArgStructure(ac, av);
            getArgVals(ac, av);
         }
    }
    
    // then read any overrides
    cout << "Reading command line parameters\n";
    getArgVals(count, args);
    check();
}    
