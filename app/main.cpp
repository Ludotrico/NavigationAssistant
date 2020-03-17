// main.cpp
//
// ICS 46 Winter 2020
// Project #5: Rock and Roll Stops the Traffic
//
// This is the program's main() function, which is the entry point for your
// console user interface.

#include "TrafficAssistant.hpp"
#include "RoadMapReader.hpp"
#include "TripReader.hpp"
#include "RoadMapWriter.hpp"


int main()
{
    //Create RoadMap Digraph and Trips vector
    InputReader in(std::cin);
    RoadMapReader mReader;
    RoadMap roadMap = mReader.readRoadMap(in);
    TripReader tripReader;
    std::vector<Trip> trips = tripReader.readTrips(in);

    //Display roadMap
    RoadMapWriter mapWriter;
    mapWriter.writeRoadMap(std::cout, roadMap);

    //Create assistnat object
    Assistant assistant(roadMap, trips);
    //Run Traffic Assistant
    assistant.run();


    return 0;
}

