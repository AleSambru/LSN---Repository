#ifndef __functions_h__
#define __functions_h__
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"


using namespace std;


// This file contains all the function and classes used in the exercise: 
/* INDEX:
 1 - Function setup
 2 - Class Position
 3- 
*/

/*
Position class:
*/
class Position {

public:

    // constructors
    Position() {
        m_x = 0;
        m_y = 0;
        m_z = 0;
    };

    Position(double x, double y, double z) {
        m_x = x;
        m_y = y;
        m_z = z;

    };

    // distructor
    ~Position() {};

    // methods
    // return cartesians coordinates
    double getX() const { return m_x; };

    double getY() const { return m_y; };

    double getZ() const { return m_z; };

    // return spherical coordinates
    double getR() const { return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2)); };

    double getPhi() const { return atan(m_y / m_x); };

    double getTheta() const { return acos(m_z / getR()); };

    double getRho() const { sqrt(pow(m_x, 2) + pow(m_y, 2)); };     // ray of the cilindric coordinates

    // calcuate the distance between the point and another point Position
    double Distance(const Position &) const {
        return sqrt(pow(getX() - m_x, 2) + pow(getY() - m_y, 2) + pow(getZ() - m_z, 2));
    }; // distanza da un altro punto

    void discrete_step() {};

    void continuous_step() {};

private:

    double m_x, m_y, m_z;

};

#endif 




