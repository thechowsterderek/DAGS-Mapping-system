// Digraph.hpp

//Digraph implementation of a mapping system
#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <queue>
#include <iostream>



// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException
{
public:
    DigraphException(const std::string& reason): reason_{reason} { }

    std::string reason() const { return reason_; }

private:
    std::string reason_;
};



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a template
// struct.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a template struct.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
    int vertexNum;
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
	class compareFunction{ //used in Priority Queue for shortestPathFunction
	public:
		bool operator()(std::pair<int,double> pOne,std::pair<int,double> pTwo){
			if(pOne.second < pTwo.second){
				return true;
			}
			else{
				return false;
			}
		}
	};
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph();

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

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
    int vertexCount() const;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const;

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

    bool vertexInMap(std::map<int,DigraphVertex<VertexInfo, EdgeInfo>>, int vertNum)const;

private:
    // Add whatever member variables you think you need here.  One
    // possibility is a std::map where the keys are vertex numbers
    // and the values are DigraphVertex<VertexInfo, EdgeInfo> objects.

	int verticesNum;
	int edgesNum;
	std::map<int,DigraphVertex<VertexInfo, EdgeInfo>> fullMap;


    // You can also feel free to add any additional member functions
    // you'd like (public or private), so long as you don't remove or
    // change the signatures of the ones that already exist.
};



// You'll need to define the member functions of your Digraph class
// template here.

template<typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo,EdgeInfo>::Digraph(){ //initialize all private variables
	verticesNum = 0;
	edgesNum = 0;
}


template<typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo,EdgeInfo>::Digraph(const Digraph& d){


	verticesNum = d.verticesNum;
	edgesNum = d.edgesNum;
	fullMap = d.fullMap;

}

template<typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo,EdgeInfo>::~Digraph(){ 
	fullMap.clear();
}

template<typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo,EdgeInfo>& Digraph<VertexInfo,EdgeInfo>::operator=(const Digraph& d){
	if(this != &d){
		fullMap.clear();
		verticesNum = d.verticesNum;
		edgesNum = d.edgesNum;
		fullMap = d.fullMap;
	}
	return *this;
}

template<typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo,EdgeInfo>::vertices()const{
	std::vector<int> verts;
	for(auto it = fullMap.begin();it!=fullMap.end();it++){//iterate through entire map for vertex numbers
		verts.push_back(it->first);
	}
	return verts;
}

template<typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo,EdgeInfo>::edges()const{
	std::vector<std::pair<int,int>> allEdges;
	std::pair<int,int> path;
	for(auto mapIter = fullMap.begin(); mapIter != fullMap.end(); mapIter++){
		for(auto it = fullMap.at(mapIter->first).edges.begin();it!=fullMap.at(mapIter->first).edges.end();it++){
			path.first = it->fromVertex;
			path.second = it->toVertex;
			allEdges.push_back(path);
		}
	}
	return allEdges;
}

template<typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo,EdgeInfo>::edges(int vertex)const{
	std::vector<std::pair<int,int>> outgoingEdges;
	std::pair<int,int> path;
	for(auto it = fullMap.at(vertex).edges.begin();it!=fullMap.at(vertex).edges.end();it++){
		path.first = it->fromVertex;
		path.second = it->toVertex;
		outgoingEdges.push_back(path);
	}
	return outgoingEdges;
}


template<typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo,EdgeInfo>::vertexInfo(int vertex) const{
	if(vertexInMap(fullMap,vertex) == false){
		throw DigraphException("Vertex Does Not Exist On The Digraph Try Again...");
	}
	else{
		return fullMap.at(vertex).vinfo;
	}

}

template<typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo,EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const{
	if(vertexInMap(fullMap,fromVertex) == false  || vertexInMap(fullMap,toVertex) == false){
		throw DigraphException("Edge Does Not Exist On The Digraph Try Again...");
	}
	else{
		EdgeInfo edge;
		for(auto it = fullMap.at(fromVertex).edges.begin(); it!= fullMap.at(fromVertex).edges.end(); it++){
			if(it->fromVertex == fromVertex && it->toVertex == toVertex){
				edge = it->einfo;
				return edge;
			}
		}
	}
	throw DigraphException("Edge Does Not Exist On The Digraph Try Again...");
}

template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo){
	if(vertexInMap(fullMap,vertex) == true){
		throw DigraphException("Vertex Already Exist Inside the Digraph Try Again ...");
	}
	else{
		DigraphVertex<VertexInfo,EdgeInfo> vertexPt;
		vertexPt.vinfo = vinfo;
		vertexPt.vertexNum = vertex;
		fullMap.emplace(vertex,vertexPt);
		verticesNum++;
	}
}

template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo){
	if(vertexInMap(fullMap,fromVertex) == false || vertexInMap(fullMap,toVertex) == false){
		throw DigraphException("One of the vertices do not exist in the Digraph");
	}
	else{
		DigraphEdge<EdgeInfo> edge;
		edge.fromVertex = fromVertex;
		edge.toVertex = toVertex;
		edge.einfo = einfo;
		for(auto it = fullMap.at(fromVertex).edges.begin(); it!= fullMap.at(fromVertex).edges.end(); it++){
			if(it->fromVertex == fromVertex && it->toVertex == toVertex){
				throw DigraphException("Edge already exist in the Digraph try again..");
			}		
		}
		fullMap.at(fromVertex).edges.push_back(edge);
		edgesNum++;
	}
}

template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::removeVertex(int vertex){
	if(vertexInMap(fullMap,vertex) == true){	
		for(auto mapIter = fullMap.begin(); mapIter != fullMap.end(); mapIter++){ 
			for(auto it = fullMap.at(mapIter->first).edges.begin();it!=fullMap.at(mapIter->first).edges.end();it++){
				if(it->toVertex == vertex){
					fullMap.at(mapIter->first).edges.erase(it);
					edgesNum--;
				}
			}
		}

		edgesNum = edgesNum - fullMap.at(vertex).edges.size();
		fullMap.at(vertex).edges.clear();
		fullMap.erase(vertex);
		verticesNum--;
		return;
	}
	else{
		throw DigraphException("Vertex does not exist in the Digraph");
	}
}

template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::removeEdge(int fromVertex, int toVertex){
	if(vertexInMap(fullMap,fromVertex) == false || vertexInMap(fullMap,toVertex) == false){
		throw DigraphException("One of the vertices do not exist in Digraph");
	}
	else{
		for(auto it = fullMap.at(fromVertex).edges.begin(); it!= fullMap.at(fromVertex).edges.end(); it++){
			if(it->fromVertex == fromVertex && it->toVertex == toVertex){
				fullMap.at(fromVertex).edges.erase(it);
				edgesNum--;
				return;
			}		
		}
		throw DigraphException("edge do not exist in Digraph");

	}
}

template<typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo,EdgeInfo>::vertexCount() const{
	return verticesNum;
}

template<typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo,EdgeInfo>::edgeCount() const{
	return edgesNum;
}

template<typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo,EdgeInfo>::edgeCount(int vertex) const{
	return fullMap.at(vertex).edges.size();
}

template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo,EdgeInfo>::isStronglyConnected() const{
	std::vector<int> verts = vertices();
	std::vector<int> vList;
	std::map<int, bool> pathCheck; //a full check if there is path for each vertex

	for(int i = 0; i < verts.size(); i++){
		pathCheck.emplace(verts[i],false); //start with all vertices connected is false
	}

	if(verts.size() > 0){
		vList.push_back(verts[0]); //only need one vertex to check if all are connected to each other
	}
	else{
		return false; //return false if there are no vertices added to Digraph
	}

	while(vList.empty() == false){ //keeps going til all vertices are went through
		int startVert = vList.back();
		vList.pop_back();

		if(pathCheck.at(startVert) == false){
			pathCheck.at(startVert) = true;

			for(auto it = fullMap.at(startVert).edges.begin(); it!=fullMap.at(startVert).edges.end(); ++it){
				for(int k = 0; k < verts.size();k++){
					if(it->toVertex == verts[k] && pathCheck.at(verts[k]) == false){
						vList.push_back(verts[k]);
					}
				}
			}
		}
	}

	bool isConnected = false;
	for(auto it = pathCheck.begin(); it != pathCheck.end(); it++){
		if(it->second == true){
			isConnected = true; //if just one path is not connected 
		}
		else{
			isConnected = false;
			break;
		}
	}
	return isConnected;
}

template<typename VertexInfo,typename EdgeInfo>
std::map<int,int> Digraph<VertexInfo,EdgeInfo>::findShortestPaths(int startVertex, std::function<double(const EdgeInfo&)> edgeWeightFunc)const{

	std::vector<int> allVertices = vertices(); //all verts
    std::map<int,bool> shortestPathFound; //check map
    std::map<int, double> distanceMeasured;
    std::map<int, int> prevVertex;
    std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,int>>,compareFunction> priorQueue;//priorty Queue 

	if(this->isStronglyConnected() == false){
		prevVertex.emplace(startVertex, startVertex);
		for(int i = 0; i < allVertices.size(); i++){
			prevVertex.emplace(allVertices[i],allVertices[i]);
		}
		return prevVertex;
	}
    const int MAX = std::numeric_limits<int>::max();

    distanceMeasured.emplace(startVertex,0);
    prevVertex.emplace(startVertex, startVertex);

    for(int i=0;i< allVertices.size();i++){
        shortestPathFound.emplace(allVertices[i], false);//put all vertices at the end of the map in order

        if(allVertices[i] != startVertex){
            distanceMeasured.emplace(allVertices[i],MAX);
            prevVertex.emplace(allVertices[i],NULL);
        }
    }
    std::vector<std::pair<int,int>> vertexEdges;
    std::pair<int,double> mainVert;
    mainVert.first = startVertex;
    mainVert.second = 0;

    priorQueue.push(mainVert);

    while(!priorQueue.empty()){
        int v = priorQueue.top().first;
        priorQueue.pop();

        if(shortestPathFound.at(v) == false){

            shortestPathFound.at(v) = true;
            double weight;

            for(auto mapIter = fullMap.at(v).edges.begin(); mapIter!= fullMap.at(v).edges.end(); ++mapIter){
                weight = edgeWeightFunc(mapIter -> einfo);
                for(int l=0; l< allVertices.size(); l++){

                    if(mapIter-> toVertex == allVertices[l] && shortestPathFound.at(allVertices[l]) == false){

                        int currVert = allVertices[l];
                        if(distanceMeasured.at(currVert) > distanceMeasured.at(v) + weight){
                            distanceMeasured.at(currVert) = distanceMeasured.at(v) + weight;
                            prevVertex.at(currVert) = v;
                            mainVert.first = currVert;
                            mainVert.second = distanceMeasured.at(currVert);
                            priorQueue.push(mainVert);
                        }
                    }
                }
            }
        }
    }
    return prevVertex;
}


template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo,EdgeInfo>::vertexInMap(std::map<int,DigraphVertex<VertexInfo, EdgeInfo>>, int vertNum)const{
	if(fullMap.find(vertNum) == fullMap.end()){
		return false;
	}
	else{
		return true;
	}
}

#endif // DIGRAPH_HPP

