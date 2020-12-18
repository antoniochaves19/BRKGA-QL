#include "Decoder.h"

// Sort TSol by random-keys
bool sortByRk(const TVecSol &lhs, const TVecSol &rhs) { return lhs.rk < rhs.rk; }

// Sort TSol by sol
//bool sortBySol(const TVecSol &lhs, const TVecSol &rhs) { return lhs.sol < rhs.sol; }

TSol Decoder(TSol s, int n, std::vector < std::vector <double> > dist)
{
    // save the random-key sequence of current solution 
    TSol temp = s;

    int numDecoders = 5;

    // create a initial solution of the problem
    s.fo = 0;
    for (int i=0; i<n; i++)
        s.vec[i].sol = i;

    int dec = ceil(s.vec[n].rk*numDecoders);
    //printf("\n%d (%.2lf)", dec, s.vec[n].rk);

    switch (dec)
    {
        case 1: // sort decoder
            s = Dec1(s, n, dist);
            break;

        case 2: // 2-Opt decoder
            s = Dec2(s, n, dist);
            break;
        
        case 3: // Cheapest Insertion
            s = Dec3(s, n, dist);
            break;

        case 4: // k-Farthest Insertion
            s = Dec4(s, n, dist);
            break;

        case 5: // k-Nearest Insertion
            s = Dec5(s, n, dist);
            break;
        
        default:
            break;
    }

    // calculate objective function
    s.fo = 0;
    for (int i=0; i<n; i++)
    {
        s.fo += dist[s.vec[i%n].sol][s.vec[(i+1)%n].sol];
    }

    // return initial random-key sequence and maintain the solution sequence
    for (int i=0; i<n; i++)
    {
        s.vec[i].rk = temp.vec[i].rk;
    }
    
    return s;
}

TSol Dec1(TSol s, int n, std::vector < std::vector <double> > dist) // sort
{
    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);   

    return s;
}

TSol Dec2(TSol s, int n, std::vector < std::vector <double> > dist) // 2-Opt
{
    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    int t = 0, i = 0, j = 0, Mi1= 0, Mj = 0;

    float foOpt = 0;

    t = n; // use a circular list
    for (i=0; i < t; i++)
    {
        j = i + 2;
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
                aux = s.vec[p1%t].sol;
                s.vec[p1%t].sol = s.vec[p2%t].sol;
                s.vec[p2%t].sol = aux;

                p1 = (p1==t-1)?0:p1+1;
                p2 = (p2 == 0)?t-1:p2-1;
            }
        }
        j++;
        }//while
    }//for
    return s;
}

TSol Dec3(TSol s, int n, std::vector < std::vector <double> > dist) // Cheapest Insertion
{
    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    TVecSol aux = s.vec[n];

    // order list of candidates
    TSol sC = s;

    // partial route with three points
    s.vec.resize(3);

    // construct a solution with cheapest insertion
    for (int i = 3; i<n; i++)
    {
        // find the cheapest position to insert the i-th point of sC
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between i-1 and 0
                costInsertion = dist[s.vec[j].sol][sC.vec[i].sol] + dist[sC.vec[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between i and i+1
                costInsertion = dist[s.vec[j].sol][sC.vec[i].sol] + dist[sC.vec[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC.vec[i]);
    }

    // last RK
    s.vec.push_back(aux);
    return s;
}

TSol Dec4(TSol s, int n, std::vector < std::vector <double> > dist) // k-Farthest Insertion
{
    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    // order list of candidates
    std::vector <TVecSol> sC = s.vec;
    sC.erase(sC.begin()+n); // apagar ultima chave de sC
    TVecSol aux = s.vec[n]; // copiar ultima chave de rk

    //TSol temp = s;

    // partial route with one point
    s.vec.clear();
    s.vec.push_back(sC[0]);
    sC.erase(sC.begin());

    // construct a solution with k farthest insertion
    while (!sC.empty())
    {
        // find the point i farthest from the partial route into the k first points of sC
        int i = 0;
        double costFarthest = -INFINITO;

        for (int k=0; k<3 && k<sC.size(); k++)
        {
            for (int j = 0; j<s.vec.size(); j++)
            {
                if (dist[s.vec[j].sol][sC[k].sol] > costFarthest)
                {
                    costFarthest = dist[s.vec[j].sol][sC[k].sol];
                    i = k;
                }
            }
        }

        // find the cheapest position to insert the point i into the partial route
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between n-1 and 0
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between j and j+1
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC[i]); //

        // erase the i-th point of the sC list
        sC.erase(sC.begin()+i);
    }

    // last random-key
    s.vec.push_back(aux);

    return s;
}

TSol Dec5(TSol s, int n, std::vector < std::vector <double> > dist) // k-Nearest Insertion
{
    // sort random-key vector 
    sort(s.vec.begin(), s.vec.end()-1, sortByRk);

    // order list of candidates
    std::vector <TVecSol> sC = s.vec;
    sC.erase(sC.begin()+n); // apagar ultima chave de sC
    TVecSol aux = s.vec[n]; // copiar ultima chave de rk

    //TSol temp = s;

    // partial route with one point
    s.vec.clear();
    s.vec.push_back(sC[0]);
    sC.erase(sC.begin());

    // construct a solution with k farthest insertion
    while (!sC.empty())
    {
        // find the point i nearest from the partial route into the k first points of sC
        int i = 0;
        double costNearest = INFINITO;

        for (int k=0; k<3 && k<sC.size(); k++)
        {
            for (int j = 0; j<s.vec.size(); j++)
            {
                if (dist[s.vec[j].sol][sC[k].sol] < costNearest)
                {
                    costNearest = dist[s.vec[j].sol][sC[k].sol];
                    i = k;
                }
            }
        }

        // find the cheapest position to insert the point i into the partial route
        int bestPosition = 0;
        float costBest = INFINITO;
        float costInsertion = 0;
        for (int j = 0; j<s.vec.size(); j++)
        {
            if (j == s.vec.size()-1)
            {
                // cost to insert between n-1 and 0
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[0].sol] - dist[s.vec[j].sol][s.vec[0].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
            else
            {
                // cost to insert between j and j+1
                costInsertion = dist[s.vec[j].sol][sC[i].sol] + dist[sC[i].sol][s.vec[j+1].sol] - dist[s.vec[j].sol][s.vec[j+1].sol];
                if (costInsertion < costBest)
                {
                    costBest = costInsertion;
                    bestPosition = j;
                }
            }
        }

        // insert the i-th point in the cheapest position
        s.vec.insert(s.vec.begin()+bestPosition+1,sC[i]); //

        // erase the i-th point of the sC list
        sC.erase(sC.begin()+i);
    }

    // last random-key
    s.vec.push_back(aux);

    return s;
}
