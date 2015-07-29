//
// Created by Andrew Pontzen on 29/07/15.
//

#ifndef IC_CHANGECWDWHILEINSCOPE_H
#define IC_CHANGECWDWHILEINSCOPE_H

#include <string>

class ChangeCwdWhileInScope {
protected:
    std::string old;

public:
    ChangeCwdWhileInScope(std::string newFolder);

    ~ChangeCwdWhileInScope();
};

std::string getDirectoryName(std::string full);


#endif //IC_CHANGECWDWHILEINSCOPE_H
