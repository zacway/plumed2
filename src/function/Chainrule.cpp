/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2022 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Function.h"
#include "ActionRegister.h"
#include "core/ActionAtomistic.h"
#include "core/PlumedMain.h"

#include <cstdio>
#include <string>
#include <iostream>
namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION CHAINRULE
/*
Calculate the time derivative of a variable by calculating a chain rule.

The functional form of this function is
\f[
C=\sum_{i=1}^{N_{atoms}} \sum_{j=1}^{3}} \frac{\partial arg}{\partial x_j^{[i]}} \frac{ \partial x_j^{[i]}}{\partial t}
\f]

It is possible to calculate the chain rule with respect to the force.
\f[
C=\sum_{i=1}^{N_{atoms}} \sum_{j=1}^{3}} \frac{\partial arg}{\partial x_i{[j]}} \frac{ Fx_i{[j]}}{\partial t}
\f]


\par Examples

The following input tells plumed to print the derivative of the distance between atoms 1 and 2.
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dAB

ca: CHAINRULE ARG=dAB PERIODIC=NO VAR=SPEED

PRINT ARG=ca FILE=out.file 
\endplumedfile
(See also \ref PRINT and \ref DISTANCE).




*/
//+ENDPLUMEDOC


class ChainRule :
  public Function ,
  public PLMD::ActionAtomistic
{
  std::string var;
  std::vector<Vector> forces;
  std::vector<Vector> positions;
  std::vector<Vector> positionsTmp;
  //to do an action only on the first call of calculate
  bool first;
public:
  explicit ChainRule(const ActionOptions&);
  void calculate() override;
  bool checkNeedsGradients()const override {return true;};
  void lockRequests() override;
  void unlockRequests() override;
  static void registerKeywords(Keywords& keys);
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override;




};
inline
void ChainRule::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

inline
void ChainRule::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

PLUMED_REGISTER_ACTION(ChainRule,"CHAINRULE")

void ChainRule::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");keys.use("PERIODIC");
  keys.add("optional","VAR","the parameters the chainrule will derive with respect to ");
 keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

ChainRule::ChainRule(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  Function(ao),
  first(true)
{
  parse("VAR",var);
  if (var.size() ==0)
    var="FORCE";


  
  std::vector<AtomNumber> Atoms;
  // TODO PROBLEM : Gradients are not calculated
  for(unsigned i=0; i<getNumberOfArguments(); ++i)    {
    std::cout<<getPntrToArgument(i)->getName()<<std::endl;
    getPntrToArgument(i)->setGradients();
    for(const auto & p : getPntrToArgument(i)->getGradients()) {
    Atoms.push_back(p.first);
  }
  }  
      requestAtoms(Atoms);
  addValueWithDerivatives();
  checkRead();
}

void ChainRule::calculate() {

  double result=0.0;
  //allows to get informations about the atoms
  Atoms& atoms(plumed.getAtoms());    
  for(const auto & p : getPntrToArgument(0)->getGradients()) {
    AtomNumber iatom=p.first;
    if (var=="FORCE"){
      atoms.getLocalMDForces(forces);
      Vector force = forces[iatom.index()];
      //chainrule for each atom force
      for (int i=0;i<3;i++) result+= p.second[i]*force[i];          
    }
    else if (var =="SPEED") {
      //initialization, the first result is 0.00 due to Numerical derivatives
      if (first==true){
              std::cout<< "GRADIENTS : "<<       (getPntrToArgument(0)->getPntrToAction())->isOptionOn("GRADIENTS")<<std::endl;
              std::cout<< "checkNeedsGradients : "<<       (getPntrToArgument(0)->getPntrToAction())->checkNeedsGradients()<<std::endl;

              positions=getPositions();
              first=false;
              return;
      }

      positionsTmp=getPositions();
      //numerical derivative of the position
      Vector dposition=(positionsTmp[iatom.index()]-positions[iatom.index()])/atoms.getTimeStep();
      //chainrule for each atom speed
      std::cout << p.second[0] << std::endl;

      for (int i=0;i<3;i++) result+= p.second[i]*dposition[i];
  
      positions=positionsTmp;
    } 
}
setValue(result);
}

void ChainRule::calculateNumericalDerivatives( ActionWithValue* a ) {
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
    }
  }
}

}
}




