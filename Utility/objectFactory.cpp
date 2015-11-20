#include "objectFactory.h"
#include "myexcept.h"
#include "../Problem/problem.h"
#include "../Algorithm/Algorithm.h"
classFactory::constructorProblem classFactory::constructProblem(const string &s){
		ClassMapProblem::iterator i=m_theMapProblem.find(s);
		if(i==m_theMapProblem.end()) throw myException("class not exists@constructorProblem::constructProblem");
		return i->second;
}

classFactory::constructorAlgorithm classFactory::constructAlgorithm(const string &s){
	ClassMapAlgorithm::iterator i=m_theMapAlgorithm.find(s);
		
	if(i==m_theMapAlgorithm.end()) throw myException("class not exists@constructorAlgorithm::constructAlgorithm");
	return i->second;
}