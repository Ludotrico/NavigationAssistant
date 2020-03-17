// Digraph.hpp
//
// ICS 46 Winter 2020
// Project #5: Rock and Roll Stops the Traffic
//
// This header file declares a class template called Digraph, which is
// intended to implement a generic directed graph.  The implementation
// uses the adjacency lists technique, so each vertex stores a linked
// list of its outgoing edges.
//
// Along with the Digraph class template is a class DigraphException
// and a couple of utility structs that aren't generally useful outside
// of this header file.
//
// In general, directed graphs are all the same, except in the sense
// that they store different kinds of information about each vertex and
// about each edge; these two types are the type parameters to the
// Digraph class template.

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include "MinHeap.hpp"
#include "Dijkstra.hpp"
#include <exception>
#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <limits>
#include <string>



// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException : public std::runtime_error
{
public:
    DigraphException(const std::string& reason);
};


inline DigraphException::DigraphException(const std::string& reason)
    : std::runtime_error{reason}
{
}



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a struct template.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a struct template.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};



// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d) noexcept;

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph() noexcept;

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d) noexcept;

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const noexcept;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const noexcept;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the precedessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;

    //Calculates Dijkstra's algorithm on a vertex and returns a vector of Dijkstra objects
    std::vector<std::vector<Dijkstra> >* getShortestPathVector(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc, std::vector<std::vector<Dijkstra> >* shortestPaths) const;


private:
    std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >* allVertices;

    //Destructor helper
    void deleteGraph();
    //Returns true if vertex exists
    bool isValidVertex(int vertex) const;
    //Returns linked list index of an edge
    int findEdge(int fromV, int toV);
    //isStrongConnected helper that returns count of how many vertices visited
    int depthFirstTraversal(int startV, std::map<int, bool>& unvisited, int& count) const;
    //Returns second-dimensional index for a vertex in shortestPaths
    int findIndex(int firstD, int toVertex, std::vector<std::vector<Dijkstra> >* shortestPaths) const;

    //findShortestPaths helper that returns distance sum from a vertex to another
    double findDistanceSum(int firstD, int fromV, int toV, std::function<double(const EdgeInfo&)> edgeWeightFunc, std::vector<std::vector<Dijkstra> >* shortestPaths) const;
    
    //Performs Dijkstra's Algorithm
    void dijkstraAlgorithm(
        int& startVertex,
        std::function<double(const EdgeInfo&)>& edgeWeightFunc, std::vector<std::vector<Dijkstra> >* shortestPaths, std::map<int, int>& result,
        std::vector<int>& unvisited, typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator& Mitr, int& index) const;
    

    

    //Stores total vertex and edge count
    int vertexAmount;
    int edgeAmount;

};


//Constructor
template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph()
{
    allVertices = new std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >;

    vertexAmount = 0;
    edgeAmount = 0;
}

//Copy constructor
template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
{
    allVertices = new std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >;
    
    *allVertices = *(d.allVertices);
    vertexAmount = d.vertexAmount;
    edgeAmount = d.edgeAmount;
}

//Copy constructor from expiring Digraph
template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d) noexcept
{
    allVertices = d.allVertices;
    vertexAmount = d.vertexAmount;
    edgeAmount = d.edgeAmount;
    d.allVertices = nullptr;
}

//Cleans up Digraph object
template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::deleteGraph()
{
    if(allVertices != nullptr) {
        delete allVertices;
    }

    
}

//Destructor
template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::~Digraph() noexcept
{
    deleteGraph();
}

//Overload = operator
template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
{
    deleteGraph();

    allVertices = new std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >;

    *allVertices = *(d.allVertices);
    vertexAmount = d.vertexAmount;
    edgeAmount = d.edgeAmount;
    
    return *this;
}

//Overload = operator with expiring Digraph
template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d) noexcept
{
    deleteGraph();

    allVertices = d.allVertices;
    vertexAmount = d.vertexAmount;
    edgeAmount = d.edgeAmount;
    d.allVertices = nullptr;

    return *this;
}

//returns a vector containing the vertex numbers of all vertices
template <typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
{
    //create map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();

    std::vector<int> v;
    //Walk through allVertices
    while (Mitr != allVertices->end()) {
        v.push_back(Mitr->first);
        Mitr++;
    }

    return v;
}

//Returns a vector containg vertex pairs of all edges
template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
{
    //create map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();

    std::vector<std::pair<int, int>> v;

    //Goes through all vertexes
    while (Mitr != allVertices->end()) {
        //create list iterator
        typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = Mitr->second.edges.begin();
        //goes through all edges
        while(Litr != Mitr->second.edges.end()) {
            v.push_back(std::make_pair(Litr->fromVertex, Litr->toVertex));
            Litr++;
        }
        
        Mitr++;
    }

    return v;
}

//Returns true if a vertex exists
template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isValidVertex(int vertex) const
{
    //create map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();

    //Walk through allVertices
    while(Mitr != allVertices->end()) {
        if(Mitr->first == vertex)
            return true;
        Mitr++;
    }

    return false;
}

//Returns a vector containing vertex pairs of all edges in given vertex
template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
{
    std::vector<std::pair<int, int>> v;
    //Check if vertex is valid
    if(!isValidVertex(vertex))
        throw DigraphException("Vertex " + std::to_string(vertex) + " does not exist.");

    //List iterator
    typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[vertex].edges.begin();

    //Go through all edges
    while(Litr != (*allVertices)[vertex].edges.end()) {
        v.push_back(std::make_pair(Litr->fromVertex, Litr->toVertex));
        Litr++;
    }

    return v;
}

//Returns vertexInfo of a given vertex
template <typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
{
    if(!isValidVertex(vertex))
        throw DigraphException("Vertex " + std::to_string(vertex) + " does not exist.");

    return (*allVertices)[vertex].vinfo;
}

//Returns edge info of a given edge
template <typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
{
    //Check if both vertices are valid
    if(!isValidVertex(fromVertex))
        throw DigraphException("fromVertex "  + std::to_string(fromVertex) + " does not exist.");
    if(!isValidVertex(toVertex))
        throw DigraphException("toVertex "  + std::to_string(toVertex) + " does not exist.");

    //list iterator
    typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[fromVertex].edges.begin();
    //Go through all edges outgoing fromVertex
    while(Litr != (*allVertices)[fromVertex].edges.end()) {
        if(Litr->toVertex == toVertex)
            return Litr->einfo;
        Litr++;
    }

    //Case that edge was not found
    throw DigraphException("Edge from Vertex " + std::to_string(fromVertex) + " to Vertex " + std::to_string(toVertex) + " does not exist.");
}

//Adds a vertex to graph
template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
{
    //Check if vertex already exists
    if((*allVertices).count(vertex) > 0)
        throw DigraphException("Vertex " + std::to_string(vertex) + " already exists.");
    
    DigraphVertex<VertexInfo, EdgeInfo> v;
    v.vinfo = vinfo;

    (*allVertices)[vertex] = v;

    vertexAmount++;
}

//Adds an edge to Digraph
template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
{
    //Check if vertices are valid
    if(!isValidVertex(fromVertex))
        throw DigraphException("fromVertex " + std::to_string(fromVertex)+ " does not exist.");
    if(!isValidVertex(toVertex))
        throw DigraphException("toVertex "  + std::to_string(toVertex) + " does not exist.");
    //Check if edge already exists
    if(findEdge(fromVertex, toVertex) != -1)
        throw DigraphException("Edge from Vertex " + std::to_string(fromVertex) + " to Vertex " + std::to_string(toVertex) + " already exists.");
    
    DigraphEdge<EdgeInfo> e;
    e.einfo = einfo;
    e.fromVertex = fromVertex;
    e.toVertex = toVertex;

    (*allVertices)[fromVertex].edges.push_front(e);

    edgeAmount++;
}

//Returns index of a given edge
template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::findEdge(int fromVertex, int toVertex)
{
    //Create list iterator
    typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[fromVertex].edges.begin();

    int index = 0;
    //Go through all edges
    while(Litr != (*allVertices)[fromVertex].edges.end()) {
        if((Litr->fromVertex == fromVertex) && (Litr->toVertex == toVertex) )
            return index;
        Litr++;
        index++;
    }

    //Never found edge
    return -1;
}

//Removes a vertex from Digraph
template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
{
    //Check if vertex does not exist
    if(!isValidVertex(vertex))
        throw DigraphException("Vertex " + std::to_string(vertex) + " does not exist.");
        
    //Delete outgoing edges
    edgeAmount -= (*allVertices)[vertex].edges.size();
    (*allVertices)[vertex].edges.clear();
    allVertices->erase(vertex);

    vertexAmount--;

    //Map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();
    //Delete incoming edges
    while(Mitr != allVertices->end()) {
        typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = Mitr->second.edges.begin();
        int index = findEdge(Mitr->first, vertex);
        if(index != -1) {
            std::advance(Litr, index);
            Mitr->second.edges.erase(Litr);
            edgeAmount--;
        }
        Mitr++;
    }

}

//Removes given edge from Digraph
template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
{
    //Check if vertices are valid
    if(!isValidVertex(fromVertex))
        throw DigraphException("fromVertex " + std::to_string(fromVertex) + " does not exist.");
    if(!isValidVertex(toVertex))
        throw DigraphException("toVertex " + std::to_string(toVertex) + " does not exist.");
    //Check if edge exists
    int index = findEdge(fromVertex, toVertex);
    if(index == -1)
        throw DigraphException("Edge from Vertex " + std::to_string(fromVertex)+ " to Vertex " + std::to_string(toVertex) + " does not exist.");
    
    //list iterator
    typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[fromVertex].edges.begin();
    //Move iterator forward to edge position
    std::advance(Litr, index);

    (*allVertices)[fromVertex].edges.erase(Litr);
    edgeAmount--;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::vertexCount() const noexcept
{
    return vertexAmount;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount() const noexcept
{
    return edgeAmount;
}

//Returns edgecount of a given vertex
template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
{
    //Check if vertex is valid
    if(!isValidVertex(vertex))
        throw DigraphException("Vertex " + std::to_string(vertex) + " does not exist.");
    
    return (*allVertices)[vertex].edges.size();
}

//DFT returns count of vertices visited
template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::depthFirstTraversal(int startV, std::map<int, bool>& visited, int& count) const
{
    //Mark startV as visited
    visited[startV] = true;
    count++;
    
    //Visit all of startV's connected, unvisited vertices
    typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[startV].edges.begin();
    while(Litr != (*allVertices)[startV].edges.end()) {
        if( visited[Litr->toVertex] == false )
            depthFirstTraversal(Litr->toVertex, visited, count);
        Litr++;
    }

    return count;
}

//Returns true if graph is strongly connected
template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
{
    //Stores visited vertices
     std::map<int, bool> visited;
    
    //Map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();
    
    //Depth first traverse every vertex
    while(Mitr != allVertices->end()) {
        //map iterator
        typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr2 = allVertices->begin();
        //Fill visited map with all vertices, false because they have not been visited
        while(Mitr2 != allVertices->end()) {
            visited[Mitr2->first] = false;
            Mitr2++;
        }
        //DFT
        int count = 0;
        if(depthFirstTraversal(Mitr->first, visited, count) != vertexAmount) {
            std::cout << Mitr->first << " is disconnected!" << std::endl;
            return false;
        }
        Mitr++;
    }

    return true;
}

//Returns index (second-dim) of vertex in shortestPaths
template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::findIndex(int firstD, int vertex, std::vector<std::vector<Dijkstra> >* shortestPaths) const
{
    for(int i = 0; i < (*shortestPaths)[firstD].size(); i++) {
        if((*shortestPaths)[firstD][i].v == vertex)
            return i;
    }
    return -1;
}

//Returns distance sum from all vertices to get to toVertex
template <typename VertexInfo, typename EdgeInfo>
double Digraph<VertexInfo, EdgeInfo>::findDistanceSum(int firstD, int fromVertex, int toVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc, std::vector<std::vector<Dijkstra> >* shortestPaths) const
{
    //Set totalDist
    int secondD = findIndex(firstD, fromVertex, shortestPaths);
    double totalDist = (*shortestPaths)[firstD][secondD].d;
    secondD = findIndex(firstD, toVertex, shortestPaths);
    
    //list iterator
    typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[fromVertex].edges.begin();
    //Go through all edges
    while(Litr != (*allVertices)[fromVertex].edges.end()) {
        if(Litr->toVertex == toVertex) {
            break;
        }
        Litr++;
    }

    //Update totalDist
    totalDist += edgeWeightFunc(Litr->einfo);

    return totalDist;
}


//Performs Dijkstra's Algorithm
template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::dijkstraAlgorithm(
    int& startVertex,
    std::function<double(const EdgeInfo&)>& edgeWeightFunc, std::vector<std::vector<Dijkstra> >* shortestPaths, std::map<int, int>& result,
    std::vector<int>& unvisited, typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator& Mitr, int& index) const
{
    //Create min heap
    MinHeap h(vertexAmount);
    //Prepare shortestPaths for push_backs
    shortestPaths->push_back(std::vector<Dijkstra>());

    //Initialize shortestPaths vector of Dijkstra objects and intialize min heap
    while(Mitr != allVertices->end()) {
        if(Mitr->first != startVertex) {
            (*shortestPaths)[index].push_back(Dijkstra(Mitr->first, std::numeric_limits<int>::max(), -1));
            h.insertKey(Dijkstra(Mitr->first, std::numeric_limits<int>::max(), -1));
        } else {
            (*shortestPaths)[index].push_back(Dijkstra(startVertex, 0, startVertex));
            h.insertKey(Dijkstra(startVertex, 0, startVertex));
        }
        Mitr++;
    }
    
    
    //====== Dijkstra's Loop ========

    //While vertices remain unvisited
    while(unvisited.size() > 0) {
        //Get Dijkstra with the smallest distance
        Dijkstra curr = h.extractMin();
        //list iterator
        typename std::list<DigraphEdge<EdgeInfo> >::iterator Litr = (*allVertices)[curr.v].edges.begin();

        //Goes to all of curr's connected, unvisited vertices
        while(Litr != (*allVertices)[curr.v].edges.end()) {
            //case that toVertex is unvisited
            if(  std::count(unvisited.begin(), unvisited.end(), Litr->toVertex) > 0 ) {
                //Calculate distance
                double distSum = findDistanceSum(index, curr.v, Litr->toVertex, edgeWeightFunc, shortestPaths);
                //find index in shortestPaths
                int vIndex = findIndex(index, Litr->toVertex, shortestPaths);
                //case that calculated distSum is less than known distance
                if( distSum < (*shortestPaths)[index][vIndex].d) {
                    //Update shortestPaths cache
                    (*shortestPaths)[index][vIndex].d = distSum;
                    (*shortestPaths)[index][vIndex].pV = curr.v;
                    //Update min heap
                    int harrIndex = h.findElement(Litr->toVertex);
                    //Decrease priority value
                    h.decreaseKey(harrIndex, distSum);
                    harrIndex = h.findElement(Litr->toVertex);
                    //change previous vertex
                    h.changePrevVertex(harrIndex, curr.v);
                }
            }
            Litr++;
        }
        //Delete curr vertex from unvisited vector
        for (std::vector<int>::iterator i = unvisited.begin(); i != unvisited.end(); ++i) {
            if (*i == curr.v ) {
                unvisited.erase(i);
                break;
            }
        }
    }

    //Go through shortestPaths vector and convert to map for output
    for(int i = 0; i < (*shortestPaths)[index].size(); i++ )
        result[(*shortestPaths)[index][i].v] = (*shortestPaths)[index][i].pV;
    
}


//Returns map of Dijkstra's algorithm
template <typename VertexInfo, typename EdgeInfo>
std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc) const
{
    std::map<int, int> result;
    
    std::vector<std::vector<Dijkstra> >* shortestPaths = new std::vector<std::vector<Dijkstra> >;
    
    //Vector stores unvisited vertices
    std::vector<int> unvisited;


    //map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();
    //Set all vertices as unvisited
    while(Mitr != allVertices->end()) {
        unvisited.push_back(Mitr->first);
        Mitr++;
    }

    Mitr = allVertices->begin();
    int index = 0;
    
    dijkstraAlgorithm(startVertex, edgeWeightFunc, shortestPaths, result, unvisited, Mitr, index);

    shortestPaths->clear();
    delete shortestPaths;

    return result;
}


//Returns map of Dijkstra's algorithm
template <typename VertexInfo, typename EdgeInfo>
std::vector<std::vector<Dijkstra> >* Digraph<VertexInfo, EdgeInfo>::getShortestPathVector(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc, std::vector<std::vector<Dijkstra> >* shortestPaths) const
{
    std::map<int, int> result;

    //Vector stores unvisited vertices
    std::vector<int> unvisited;


    //map iterator
    typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo> >::iterator Mitr = allVertices->begin();
    //Set all vertices as unvisited
    while(Mitr != allVertices->end()) {
        unvisited.push_back(Mitr->first);
        Mitr++;
    }

    //Initialize shortestPaths vector of Dijkstra structs and min heap
    Mitr = allVertices->begin();
    int index = shortestPaths->size();
    
    dijkstraAlgorithm(startVertex, edgeWeightFunc, shortestPaths, result, unvisited, Mitr, index);
    
    return shortestPaths;
}

#endif // DIGRAPH_HPP



