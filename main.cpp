// main.cpp

// This is the program's main() function, which is the entry point for interface


#include "InputReader.hpp"
#include "RoadMap.hpp"
#include "TripReader.hpp"
#include <iostream>
#include <vector>
#include "RoadMapReader.hpp"
#include <algorithm>
#include <math.h>


double miles(RoadSegment i){ //return miles for edge
	return i.miles;
}
double times(RoadSegment i){//mph for edge
	return i.milesPerHour;
}
std::string timeRead(double b){//convert the time to a human readable format.
	if(b>1){
		int hours = b;
		double mins = fmod(b*60,60);
		double seconds = fmod((mins*60),60);
		int realmin = mins;
		int reals = seconds;
		std::string total = std::to_string(hours)+  std::to_string(realmin) + " Minutes"+" " +std::to_string(reals)+"Seconds";
		return total;
	}
	else{
		double mins = b*60;
		double seconds = fmod(mins*60,60);
		int realm = mins;
		int reals = seconds;
		std::string total = std::to_string(realm) + " Minutes"+" " +std::to_string(reals) +" Seconds";
		return total;
	}
}

int main()
{	
	std::string temp;
	RoadMap map;
	std::vector<Trip> tripList;//list of all trips
	std::map<int,int> allPathways;
	std::vector<int> entireTrip; //The entire  pathway of trip


	InputReader d(std::cin);
	RoadMapReader r;
	map = r.readRoadMap(d);
	if(map.isStronglyConnected() == true){//check if its strongly connected
		TripReader t;
		tripList = t.readTrips(d);
		for(int k =0; k< tripList.size();k++){
			if(tripList[k].metric == TripMetric::Distance){
				double total = 0.0;
				std::function<double(RoadSegment)> f = miles;
				allPathways = map.findShortestPaths(tripList[k].startVertex,f);//pathway of entire trip
				while(allPathways.find(tripList[k].endVertex)->second != tripList[k].startVertex){ //follow the list back to the start
					entireTrip.push_back(tripList[k].endVertex);//add the next step to the list
					tripList[k].endVertex = allPathways.find(tripList[k].endVertex)->second; //move to the next step
				}
				entireTrip.push_back(tripList[k].startVertex);//add the starting vertex to the list of steps
				std::reverse(entireTrip.begin(),entireTrip.end());//reverse the list, so starting vertex is first
				for (int i = 0;i<entireTrip.size();i++){//readable format
					if (i==0){
					std::cout <<"Begin at " << map.vertexInfo(entireTrip[i]) << std::endl;
				}
				else{
					std::cout <<"Continue to " << map.vertexInfo(entireTrip[i]) <<"("<< map.edgeInfo(entireTrip[i],entireTrip[i+1]).miles<<" miles)" << std::endl;
					total += map.edgeInfo(entireTrip[i],entireTrip[i+1]).miles;
				}

			}
			std::cout<< "Total Miles = " << total << " miles" << std::endl;
			std::cout<< std::endl;

		}
		else{
			double totalTime = 0;

			std::function<double(RoadSegment)> f = times;//same thing as above but this time calculate the time it takes
			allPathways = map.findShortestPaths(tripList[k].startVertex,f);//then reorganize the time to readable format
			while(allPathways.find(tripList[k].endVertex)->second != tripList[k].startVertex){
				entireTrip.push_back(tripList[k].endVertex);
				tripList[k].endVertex = allPathways.find(tripList[k].endVertex)->second; 
			}
			entireTrip.push_back(tripList[k].startVertex);
			std::reverse(entireTrip.begin(),entireTrip.end());
			for (int i = 0;i<entireTrip.size();i++){
				if (i==0){
					std::cout <<"Begin at " << map.vertexInfo(entireTrip[i]) << std::endl;
				}
				else{
					double temp;
					temp = map.edgeInfo(entireTrip[i],entireTrip[i+1]).miles/map.edgeInfo(entireTrip[i],entireTrip[i+1]).milesPerHour;
					std::cout <<"Continue to " << map.vertexInfo(entireTrip[i]) <<"("<< map.edgeInfo(entireTrip[i],entireTrip[i+1]).miles<< "@"<<
						map.edgeInfo(entireTrip[i],entireTrip[i+1]).milesPerHour<<"mph"<< " "<< timeRead(temp) << ")" << std::endl;
					totalTime += temp;
				
				}

			}
			std::cout<< "Total Time = " << timeRead(totalTime) << std::endl;
			std::cout<< std::endl;

			}
		}
	}
	else{
		std::cout << "Disconnected Map, not all locations reachable to each other" <<std::endl; //output disconnected in case
	}
}	


