#ifndef TREEEVENTMAP_H
#define TREEEVENTMAP_H
//#include "treeeventnode.h"
#include <map>
#include <R.h>
#include <Rinternals.h>
#include "eventseq.h"
//using namespace std;

//Forward declaration
class TreeEventNode;
//Iterateur pour la map
typedef std::map<int,TreeEventNode*>::iterator TreeEventMapIterator;
//D�finition de type pour aleger le code
//Une map sur le type d'�v�nements, et la classe �v�nement
class TreeEventMap: public std::map<int,TreeEventNode*> {

public:
TreeEventMap():std::map<int,TreeEventNode*>() {}
    ~TreeEventMap() {}
    void simplifyTreeMap(const int &minSup);
    void print(const int & prof, const bool& isbrother);
    int countSubsequence(int minSup);
    void getSubsequences(SEXP result,int * support, Sequence *s2, int *index,const double &step, SEXP classname,EventDictionary * ed);
    void clearAllPointers();
    void clearSupport();

};

#endif // TREEEVENTMAP_H
