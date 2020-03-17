#include "Digraph.hpp"


Dijkstra::Dijkstra(int vertex, double distance, int previousVertex, bool b) 
{
    v = vertex;
    d = distance;
    pV = previousVertex;

    isDist = b;
};

Dijkstra::Dijkstra() : v{-1}, d{-1}, pV{-1}, isDist{true} 
{
}

//Operator overloading for priority queue
bool operator<(const Dijkstra& d1, const Dijkstra& d2) 
{ 
    //returns true if d1 is a smaller distance away
    return d1.d < d2.d; 
} 

//Operator overloading for priority queue
bool operator>(const Dijkstra& d1, const Dijkstra& d2) 
{ 
    //returns true if d1 is a bigger distance away
    return d1.d > d2.d; 
} 
