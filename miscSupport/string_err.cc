
#include <iostream.h>
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




