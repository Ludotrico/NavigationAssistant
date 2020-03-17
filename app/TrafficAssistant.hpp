#ifndef TRAFFICASSISTANT_HPP
#define TRAFFICASSISTANT_HPP

#include "RoadMap.hpp"
#include "RoadSegment.hpp"
#include "Trip.hpp"
#include "TripMetric.hpp"
#include <string>
#include <iostream>

class Assistant {
public:
    //Constructor
    Assistant(RoadMap, std::vector<Trip>);
    //Destructor
    ~Assistant();

    //Starts Traffic Assistant
    void run();

    //Processes Trip
    void processTrip(Trip t);
    
    //Prints Trip
    void printTrip(std::vector<Dijkstra> d, Trip t);

    //Processes time given an hour measurement
    std::string processTime(double t);

    //Helper to processTime
    std::string to_string(double d);


private:
    //Digraph
    RoadMap roadMap;
    
    //Vector stores all trips
    std::vector<Trip> trips;
    
    //Returns a vector of Dijktra objects, ordered sequentially from one vertex to destination
    std::vector<Dijkstra> getTrip(int fromV, int toV);
    
    //Returns index of startVertex of a given shortestPath
    int findStart(int firstD);
    
    //Finds second-dim index of vertex in shortestPaths
    int findIndex(int firstD, int vertex);
    
    //getTrip helper
    void getTripRecursively(int fromV, int toV, int firstD, std::vector<Dijkstra>& trip);
    
    //Cache for Dijkstra's algorithm
    std::vector<std::vector<Dijkstra> >* shortestPaths;
    
    //True if trip is optimizing distance, false if optimizing time
    bool dist;
};



//Constructor
Assistant::Assistant(RoadMap r, std::vector<Trip> t )
    : roadMap{r}, trips{t}
{
    shortestPaths = new std::vector<std::vector<Dijkstra> >;
}

//Destructor
Assistant::~Assistant() {
    if(shortestPaths != nullptr) {
        (*shortestPaths).clear();
        delete shortestPaths;
    }
}

//Processes each trip
void Assistant::run() {
    for(int i = 0; i < trips.size(); i++) {
        processTrip(trips[i]);
    }
}


//Helper to getTrip
int Assistant::findStart(int firstD)
{
    for(int i = 0; i < (*shortestPaths)[firstD].size(); i++) {
        if( (*shortestPaths)[firstD][i].d == 0)
            return i;
    }
    return -1;
    
}

//Returns a vector of Dijktra objects, ordered sequentially from one vertex to destination
std::vector<Dijkstra>  Assistant::getTrip(int fromV, int toV)
{
    int firstD = 0;
    //Finds trip in shortestPaths
    for(int i = 0; i < shortestPaths->size(); i++) {
        if( ((*shortestPaths)[i][findStart(i)].v == fromV) && ((*shortestPaths)[i][0].isDist != dist) ) {
            firstD = i;
            break;
        }
    }
    
    std::vector<Dijkstra> trip;
    getTripRecursively(fromV, toV, firstD, trip);

    return trip;
}

//Returns index (second-dim) of vertex in shortestPaths
int Assistant::findIndex(int firstD, int vertex)
{
    for(int i = 0; i < (*shortestPaths)[firstD].size(); i++) {
        if((*shortestPaths)[firstD][i].v == vertex)
            return i;
    }
    return -1;
}


//Returns a vector of Dijktra objects, ordered sequentially from one vertex to destination
void  Assistant::getTripRecursively(int fromV, int toV, int firstD, std::vector<Dijkstra>& trip)
{
    //Find second-dim
    int index = findIndex(firstD, toV);
    
    
    if((*shortestPaths)[firstD][index].v != fromV) {
        int pV = (*shortestPaths)[firstD][index].pV;
        //Each frame creates new frame with what vertex it came from
        getTripRecursively(fromV, pV, firstD, trip);
    }
    
    //Base case where recursion has tracked back to the start vertex
    trip.push_back((*shortestPaths)[firstD][index]);
}



//Processes a trip
void Assistant::processTrip(Trip trip) {
    dist = !(bool)trip.metric;
    bool inCache = false;
    
    //Check cache to see if startVertex has already been computed
    if((*shortestPaths).size() > 0) {
        for(int i = 0; i < shortestPaths->size(); i++) {
            if( ((*shortestPaths)[i][findStart(i)].v == trip.startVertex) &&  ((*shortestPaths)[i][0].isDist != dist) ) {
                printTrip( getTrip(trip.startVertex, trip.endVertex), trip);
                inCache = true;;
            }
        }
    }
    //Case that startVertex is not in cache
    if(!inCache) {
        //Calculate shortest path (edge weight Distance)
        if (trip.metric == TripMetric::Distance) {
            roadMap.getShortestPathVector(trip.startVertex, [](RoadSegment e)
                {
                    return e.miles;
                }, shortestPaths);
        }
        //Calculate shortest path (edge weight time)
        else {
            roadMap.getShortestPathVector(trip.startVertex, [](RoadSegment e)
                {
                    return ((double)e.miles)/e.milesPerHour;
                }, shortestPaths);
        }
        
        //Update shortestPaths to reflect what is being optimized (distance/time)
        (*shortestPaths)[(*shortestPaths).size()-1][0].isDist = !dist;
        
        //Prints trip
        printTrip( getTrip(trip.startVertex, trip.endVertex), trip);
    }
}

//Helper to processtime, fixes trailing zeros
std::string Assistant::to_string(double d) {
    std::string str = "";
    str = std::to_string(d);
    str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
    return str;
    
}



//Processes time
std::string Assistant::processTime(double time) {
    int hours = 0;
    int minutes = 0;
    double seconds = 0;
    std::string s = "";

    //Case to only use seconds
    if(time  < ((double)1/60) ) {
        seconds = time * 60 * 60;
        s += to_string((int)seconds) + "." + to_string((int) ((seconds-(int)seconds) * 10)) + ((seconds > 1)?" secs":" sec");
    }
    //Case to use minutes and seconds
    else if (time >= ((double)1/60) && time < 1) {
        minutes = time*60;
        seconds = ((time*60) - minutes) * 60;
        s += to_string(minutes) + ((minutes > 1)?" mins ":" min ") + to_string((int)seconds) + "." + to_string((int) ((seconds-(int)seconds) * 10)) + ((seconds > 1)?" secs":" sec");
    }
    //Case to use hours, minutes, and seconds
    else {
        hours = time;
        minutes = (time - hours) * 60;
        seconds = (((time - hours) * 60) - minutes) * 60;
        s += to_string(hours) + ((hours > 1)?" hrs ":" hr ") + to_string(minutes) + ((minutes > 1)?" mins ":" min ") + to_string((int)seconds) + "." + to_string((int) ((seconds-(int)seconds) * 10)) + ((seconds > 1)?" secs":" sec");
    }

    return s;
}

void Assistant::printTrip(std::vector<Dijkstra> d, Trip trip) {
    std::string s = "";
    
    //Case that trip is optimizing for distance
    if(trip.metric == TripMetric::Distance) {
        s += "Shortest distance from " + roadMap.vertexInfo(trip.startVertex) + " to " + roadMap.vertexInfo(trip.endVertex) + "\n";
        s += "  Begin at " + roadMap.vertexInfo(trip.startVertex) + "\n";
        d.erase(d.begin());
        
        //Go through each step of the trip
        for(int i = 0; i < d.size(); i++) {
            double distance = roadMap.edgeInfo(d[i].pV, d[i].v).miles;
            s += "  Continue to " + roadMap.vertexInfo(d[i].v) + " (" + to_string(distance) + ((distance > 1)?" miles)":" mile)") + "\n";
        }
        s += "Total distance: " + to_string(d[d.size()-1].d) + ((d[d.size()-1].d > 1)?" miles":" mile") + "\n\n";
    }
    //Case that trip is optimizing for time
    else {
        s += "Shortest driving time from " + roadMap.vertexInfo(trip.startVertex) + " to " + roadMap.vertexInfo(trip.endVertex) + "\n";
        s += "  Begin at " + roadMap.vertexInfo(trip.startVertex) + "\n";
        d.erase(d.begin());
        
        //Go through each step of the trip
        for(int i = 0; i < d.size(); i++) {
            double mph = roadMap.edgeInfo(d[i].pV, d[i].v).milesPerHour;
            double distance = roadMap.edgeInfo(d[i].pV, d[i].v).miles;
            s += "  Continue to " + roadMap.vertexInfo(d[i].v) + " (" + to_string(distance) + " miles @ " + to_string(mph) + "mph = " + processTime(distance/mph) + ")\n";
        }

        s += "Total time: " + processTime(d[d.size()-1].d) + "\n\n";
    }
    
    //Print trip
    std::cout << s;
}





#endif // ROADMAP_HPP
