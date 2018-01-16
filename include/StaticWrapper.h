//
// Created by Wulfix on 02/12/2017.
//

#ifndef DOCIPROJECT_INTEGRALCALCULATOR_H
#define DOCIPROJECT_INTEGRALCALCULATOR_H


class StaticWrapper {
public:
    virtual double calculateOverlap(unsigned long site1, unsigned long site2)=0;
    virtual double calculateOverlap(unsigned long site1, unsigned long site2, unsigned long site3, unsigned long site4)=0;
    virtual unsigned long getN_bf()=0;
    virtual unsigned long getN_electrons()=0;
};


#endif //DOCIPROJECT_INTEGRALCALCULATOR_H
