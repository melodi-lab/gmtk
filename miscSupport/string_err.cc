
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
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

//---------------------------------------------------------------------------


const int insert_cost=1;
const int delete_cost=1;
const int sub_cost=1;

struct from_rec
{
     int row, column, cost;
     from_rec(int r, int c, int _cost) {row=r; column=c; cost=_cost;}
     from_rec() {row=column=cost=0;}
};

//---------------------------------------------------------------------------


int score(vector<string> &reference, vector<string> &actual, bool
display=false)
{
     int max_words = max(reference.size(), actual.size());
     int cost[max_words+1][max_words+1];
     from_rec from[max_words+1][max_words+1];

     int reference_words = reference.size();
     int actual_words = actual.size();

     int i,j;
     cost[0][0] = 0;
     for (j=1; j<=reference_words; j++)
     {
          cost[0][j] = cost[0][j-1] + insert_cost;
          from[0][j] = from_rec(0,j-1,insert_cost);
     }

     for (i=1; i<=actual_words; i++)
     {
          cost[i][0] = cost[i-1][0] + delete_cost;
          from[i][0] = from_rec(i-1,0,delete_cost);
     }

     for (int i=1; i<=actual_words; i++)
          for (int j=1; j<=reference_words; j++)
          {
               int c = ((actual[i-1]==reference[j-1])?0:sub_cost);
               cost[i][j] = cost[i-1][j-1] + c;
               from[i][j] = from_rec(i-1,j-1,c);

               if (cost[i-1][j]+delete_cost < cost[i][j])
               {
                    cost[i][j] = cost[i-1][j] + delete_cost;
                    from[i][j] = from_rec(i-1,j,delete_cost);
               }

               if (cost[i][j-1]+insert_cost < cost[i][j])
               {
                    cost[i][j] = cost[i][j-1] + insert_cost;
                    from[i][j] = from_rec(i,j-1,insert_cost);
               }
          }

     if (display)
     {
          cout << "The minimum cost is: "
               << cost[actual_words][reference_words] << endl;

          const int msg_len=100;
          char msg[max_words][msg_len];
          int op = -1;
          int r=actual_words, c=reference_words;
          while (r!=0 || c!=0)
          {
                    if (from[r][c].row == r)
                    {
                              strcpy(msg[++op], "Insert ");
                              strcat(msg[op], reference[c-1].c_str());
                    }
                    else if (from[r][c].column == c)
                    {
                              strcpy(msg[++op], "Delete ");
                              strcat(msg[op], actual[r-1].c_str());
                    }
                    else
                    {
                              strcpy(msg[++op], "Match ");
                              strcat(msg[op], actual[r-1].c_str());
                              strcat(msg[op], " , ");
                              strcat(msg[op], reference[c-1].c_str());
                    }
                    strcat(msg[op], " : ");
                    char s[100];
                    sprintf(s, "%i", from[r][c].cost);
                    strcat(msg[op], s);
                    from_rec f = from[r][c];
                    r = f.row;
                    c = f.column;
          }

          while (op >= 0)
                    cout << msg[op--] << endl;
     }

     return cost[actual_words][reference_words];
}

main(int argc, char *argv[])
{
    if (argc!=3 && argc!=4) 
      {cout<<"Usage: string_err infile reference [\"show-align\"]\n"; exit(1);}
 
    /* This program takes two files that have utterances separated by #.
       It computes the string WER on an utterance by utterance basis, and
       prints out the alignments, number of errors, and overall error rate.
    */

    ifstream in(argv[1]);
    if (!in) {cerr << "Unable to open " << argv[1] << endl; exit(1);}

    ifstream ref(argv[2]);
    if (!ref) {cerr << "Unable to open " << argv[2] << endl; exit(1);}

    bool show_align = false;
    if (argc==4 && !strcmp(argv[3],"show-align"))
        show_align = true;

    int errs=0,refwrds=0;
    while (!in.eof())
    {
        string s;
        vector<string> instr, refstr;

        // read the decoded string
        while (1)
        {
            in >> s >> ws;
            if (s == "#")
               break;
            instr.push_back(s);
            assert(in.good());
        }
        assert(s=="#");

        // read the reference string
        while (1)
        {
            ref >> s >> ws;
            if (s == "#")
               break;
            refstr.push_back(s);
            assert(ref.good());
        }
        assert(s=="#");
    
        refwrds += refstr.size(); 
        errs += score(refstr, instr, show_align);    
    }
        
    cout << refwrds << " reference words\n";
    cout << errs << " errors\n";
    cout << "WER: " << float(100*errs)/refwrds << "%\n";

    in.close();
    ref.close();
}

