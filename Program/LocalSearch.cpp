#include "LocalSearch.h"

TSol LocalSearch(TSol s, int n, std::vector < std::vector <double> > dist)
{
    // we use a Random Variable Neighborhood Descent (RVND) as local search
	int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    for (int i=1; i<=4; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

    //perturbar a solucao
    if (rand(0,1) < 0.5)
    {
    }

    //printf("\nHeuristicas: ");
	while (!NSL.empty())
	{
        // current objective function
        double foCurrent = s.fo;

        // randomly choose a neighborhood
        int pos = irand(0,NSL.size()-1);
        k = NSL[pos];

        switch (k)
        {
        case 1: 
            s = LS1(s, n, dist); // 2-Opt heuristic
            break;

        case 2:
            s = LS2(s, n, dist); // node-insertion heuristic
            break;

        case 3:
            s = LS3(s, n, dist); // node-exchange heuristic
            break;

        case 4:
            s = LS4(s, n, dist); //DoubleBridge(s); // double bridge heuristic
            break;
        
        default:
            break;
        }

        // return to first neighborhood if better the current solution
        if (s.fo < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        // next neighborhood, otherwise
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
	} //end while

    // convert the solution found by the local search in a random-key vector
    /*TSol randomKey = s,
         route = s;

    // sort random-key vector 
    sort(randomKey.vec.begin(), randomKey.vec.end()-1, sortByRk); 

    for (int i=0; i<n; i++)
    {
        s.vec[route.vec[i].sol].rk = randomKey.vec[i].rk;
    }  

    s.vec[n].rk = 0.01;*/

 	return s;
}

TSol LS1(TSol s, int n, std::vector < std::vector <double> > dist) // 2-Opt
{
    int t = n, // use a circular list
        i = 0,
        j = 0,
        Mi1= 0,
        Mj = 0;

    double foOpt = 0;

    if (t > 4)
    {
        for (i=0; i < t; i++)
        {
            j = i+2;
            while (((j+1)%t) != i)
            {
                int vi  = s.vec[i].sol;
                int vi1 = s.vec[(i+1)%t].sol;
                int vj  = s.vec[j%t].sol;
                int vj1 = s.vec[(j+1)%t].sol;

                foOpt = - dist[vi][vi1]
                        - dist[vj][vj1]
                        + dist[vi][vj]
                        + dist[vi1][vj1];

                if (foOpt < 0)
                {
                    // first improvement strategy
                    Mi1 = (i+1)%t;
                    Mj  = j%t;

                    int inicio = Mi1,
                    fim = Mj;

                    int tam, p1, p2, aux;

                    if(inicio > fim)
                        tam = t - inicio + fim + 1;
                    else
                        tam = fim - inicio + 1;

                    p1=inicio;
                    p2=fim;

                    for(int k=0; k < tam/2; k++)
                    {
                        //printf("\n %d %d %d %d  =>  %d", Mi1, Mj, p1%t, p2%t, t);
                        aux = s.vec[p1%t].sol;
                        s.vec[p1%t].sol = s.vec[p2%t].sol;
                        s.vec[p2%t].sol = aux;

                        p1 = (p1==t-1)?0:p1+1;
                        p2 = (p2 == 0)?t-1:p2-1;
                    }
                    s.fo = s.fo + foOpt;
                }
                j++;
            }//while
        }//for
    }//if t > 4
    return s;
}

TSol LS2(TSol s, int n, std::vector < std::vector <double> > dist) // NodeInsertion
{
    int i = 0,
        j = 0;
    
    double foOpt;

    for (i=0; i < n; i++)
    {
        j = (i+1)%n;
        while ( ((j+1)%n) != i )
        {
            int vi  = s.vec[i].sol;
            int viP = s.vec[(i+1)%n].sol;
            int viM = 0;
            if (i == 0)
                viM = s.vec[n-1].sol;
            else
                viM = s.vec[i-1].sol;
            
            int vp  = s.vec[j%n].sol;
            int vq = s.vec[(j+1)%n].sol;

            foOpt = - dist[vp][vq]
                    - dist[viM][vi]
                    - dist[vi][viP]
                    + dist[vp][vi]
                    + dist[vi][vq]
                    + dist[viM][viP];

            if (foOpt < 0)
            {
                // first improvement strategy
                TVecSol aux;
                
                aux.sol = s.vec[i].sol;
                aux.rk = s.vec[i].rk;
                
                s.vec.insert(s.vec.begin()+((j+1)%n), aux);
                
                if (i < ((j+1)%n))
                    s.vec.erase(s.vec.begin()+i);
                else
                    s.vec.erase(s.vec.begin()+(i+1));

                s.fo = s.fo + foOpt; 
            }
            j++;
        }
    }
    return s;
}

TSol LS3(TSol s, int n, std::vector < std::vector <double> > dist) // NodeExchange
{
    int i = 0,
        j = 0;
    
    double foOpt;

    for (i=0; i < n-1; i++)
    {
        j = i+2;
        while (j < n)
        {
            if (i != 0 && j != n-1) //no exchange edge of the tour
            {
                int vi  = s.vec[i].sol;
                int viP = s.vec[i+1].sol;
                int viM = 0;
                if (i == 0)
                    viM = s.vec[n-1].sol;
                else
                    viM = s.vec[i-1].sol;
                
                int vj  = s.vec[j].sol;
                int vjM = s.vec[j-1].sol;
                int vjP = 0;
                if (j < n-1)
                    vjP = s.vec[j+1].sol;
                else
                    vjP = s.vec[0].sol;

                foOpt = - dist[viM][vi]
                        - dist[vi][viP]
                        - dist[vjM][vj]
                        - dist[vj][vjP]
                        + dist[viM][vj]
                        + dist[vj][viP]
                        + dist[vjM][vi]
                        + dist[vi][vjP];

                if (foOpt < 0)
                {
                    // first improvement strategy
                    TVecSol aux;

                    // exchange i and j
                    aux = s.vec[i];
                    s.vec[i] = s.vec[j];
                    s.vec[j] = aux;
                    
                    s.fo = s.fo + foOpt; 
                }
            }
            j++;
        }
    }
    return s;
}

TSol LS4(TSol s, int n, std::vector < std::vector <double> > dist) // OrOpt2
{
    int i = 0,
        j = 0;
    
    double foOpt;

    for (i=0; i < n-1; i++)
    {
        j = i+2;
        while ( ((j+1)%n) != i )
        {
            int vi  = s.vec[i].sol;
            int viP1 = s.vec[(i+1)%n].sol;
            int viP2 = s.vec[(i+2)%n].sol;
            int viM = 0;
            if (i == 0)
                viM = s.vec[n-1].sol;
            else
                viM = s.vec[i-1].sol;
            
            int vp  = s.vec[j%n].sol;
            int vq = s.vec[(j+1)%n].sol;

            foOpt = - dist[vp][vq]
                    - dist[viM][vi]
                    - dist[viP1][viP2]
                    + dist[vp][vi]
                    + dist[viP1][vq]
                    + dist[viM][viP2];

            if (foOpt < 0) //&& i < (j%n)
            {
                // first improvement strategy
                
                // movement
                TVecSol aux1, aux2;
                
                aux1 = s.vec[i];
                aux2 = s.vec[(i+1)%n];

                if (i < ((j+1)%n))
                {
                    if (j%n == n-1)
                    {
                        s.vec.push_back(aux1);
                        s.vec.push_back(aux2);
                    }
                    else
                    {
                        s.vec.insert(s.vec.begin()+((j+1)%n), aux2); // add vi+1
                        s.vec.insert(s.vec.begin()+((j+1)%n), aux1); // add vi
                    }

                    s.vec.erase(s.vec.begin()+((i+1)%n));
                    s.vec.erase(s.vec.begin()+i);
                }
                else if (((j+1)%n) < i)
                {
                    s.vec.erase(s.vec.begin()+((i+1)%n));   // drop vi+1
                    s.vec.erase(s.vec.begin()+(i));         // drop vi

                    s.vec.insert(s.vec.begin()+((j+1)%n), aux2); // add vi+1
                    s.vec.insert(s.vec.begin()+((j+1)%n), aux1); // add vi
                }

                //update fitness
                s.fo = s.fo + foOpt; 
            }
            j++;
        }
    }
    return s;
}

double rand(double min, double max)
{
	return ((double)(rand()%10000)/10000.0)*(max-min)+min;
    //return uniform_real_distribution<double>(min, max)(rng);
}

int irand(int min, int max)
{
	return (int)rand(0,max-min+1.0) + min;
}