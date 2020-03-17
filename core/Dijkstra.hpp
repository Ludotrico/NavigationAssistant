#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP

#include <vector>



class Dijkstra
{
public:

    //Current vertex
    int v;
    //Distance unit (could be distance, time, etc)
    double d;
    //Previous vertex
    int pV;

    //distance helper
    bool isDist;

    Dijkstra(int vertex, double distance, int previousVertex, bool b=true);
    Dijkstra();

    friend bool operator < (const Dijkstra& d1, const Dijkstra& d2);
    friend bool operator > (const Dijkstra& d1, const Dijkstra& d2);     
};



#endif // DIJKSTRA_HPP

