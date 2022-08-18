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
#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/File.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPGRADIENT
/*
Dump the gradient with respect to the all the atoms positions for one or more objects (generally CVs, functions or biases).



\par Examples

The following input instructs plumed to write a file called deriv that contains the
analytical derivatives of the distance between atoms 1 and 2.
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dAB
DUMPGRADIENT ARG=dAB STRIDE=1 FILE=deriv
\endplumedfile

(See also \ref DISTANCE)
  
The format of the output file is as follow :
\auxfile{COLVAR}
 FIELDS time    IDatom d(dAB)/dx       d(dAB)/dy       d(dAB)/dz
        0.000000 1   -0.2131822274   -0.9770124554    0.0000000000
        0.000000 2    0.2131822274    0.9770124554    0.0000000000
        1.000000 1   -0.2098034148   -0.9777435897    0.0000000000
        1.000000 2    0.2098034148    0.9777435897    0.0000000000
        2.000000 1   -0.2064099087   -0.9784656098    0.0000000000
        2.000000 2    0.2064099087    0.9784656098    0.0000000000
\endauxfile


If gradients of multiple objects are asked such as in this command :
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dAB
DISTANCE ATOMS=1,3 LABEL=dAC

DUMPGRADIENT ARG=dAB,dAC STRIDE=1 FILE=deriv
\endplumedfile

Then the output format will be as follow:

\auxfile{COLVAR}
#! FIELDS time IDatom d(dAB)/dx d(dAB)/dy d(dAB)/dz
 0.000000 1   -0.2131822274   -0.9770124554    0.0000000000
 0.000000 2    0.2131822274    0.9770124554    0.0000000000
#! FIELDS time IDatom d(dAC)/dx d(dAC)/dy d(dAC)/dz
 0.000000 1   -0.2341781290    0.9721937070    0.0000000000
 0.000000 3    0.2341781290   -0.9721937070    0.0000000000
#! FIELDS time IDatom d(dAB)/dx d(dAB)/dy d(dAB)/dz
 1.000000 1   -0.2098034148   -0.9777435897    0.0000000000
 1.000000 2    0.2098034148    0.9777435897    0.0000000000
#! FIELDS time IDatom d(dAC)/dx d(dAC)/dy d(dAC)/dz
 1.000000 1   -0.2342482310    0.9721768184    0.0000000000
 1.000000 3    0.2342482310   -0.9721768184    0.0000000000
#! FIELDS time IDatom d(dAB)/dx d(dAB)/dy d(dAB)/dz
 2.000000 1   -0.2064099087   -0.9784656098    0.0000000000
 2.000000 2    0.2064099087    0.9784656098    0.0000000000
#! FIELDS time IDatom d(dAC)/dx d(dAC)/dy d(dAC)/dz
 2.000000 1   -0.2343094385    0.9721620683    0.0000000000
 2.000000 3    0.2343094385   -0.9721620683    0.0000000000
\endauxfile



*/
//+ENDPLUMEDOC

class DumpGradient :
  public ActionPilot,
  public ActionWithArguments
{
  std::string file;
  std::string fmt;
  OFile of;
public:
  void calculate() override {}
  explicit DumpGradient(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  bool checkNeedsGradients()const override {return true;};
  void apply() override {};
  void update() override;
  ~DumpGradient();
};

PLUMED_REGISTER_ACTION(DumpGradient,"DUMPGRADIENT")

void DumpGradient::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the derivatives should be output");
  keys.add("compulsory","FILE","the name of the file on which to output the derivatives");
  keys.add("compulsory","FMT","%15.10f","the format with which the derivatives should be output");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

DumpGradient::DumpGradient(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  fmt("%15.10f")
{
  parse("FILE",file);
  if( file.length()==0 ) error("name of output file was not specified");
  parse("FMT",fmt);
  fmt=" "+fmt;
  of.link(*this);
  of.open(file);
  log.printf("  on file %s\n",file.c_str());
  log.printf("  with format %s\n",fmt.c_str());
  unsigned nargs=getNumberOfArguments();
  if( nargs==0 ) error("no arguments specified");
  (getPntrToArgument(0)->getPntrToAction())->turnOnDerivatives();
  // unsigned npar=getPntrToArgument(0)->getNumberOfDerivatives();
  // if( npar==0 ) error("one or more arguments has no derivatives");
  for(unsigned i=1; i<nargs; i++) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
  //   if( npar!=getPntrToArgument(i)->getNumberOfDerivatives() ) error("the number of derivatives must be the same in all values being dumped");
  }
  checkRead();
}


void DumpGradient::update() {
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    for(const auto & p : getPntrToArgument(i)->getGradients()) {
        AtomNumber iatom=p.first;
        of.fmtField(" %f");
        of.printField("time",getTime());
        of.printField("IDatom",(int)iatom.serial());
        of.fmtField(fmt);  
        of.printField("d("+getPntrToArgument(i)->getName()+")/dx", p.second[0]  );//
        of.printField("d("+getPntrToArgument(i)->getName()+")/dy", p.second[1]  );//
        of.printField("d("+getPntrToArgument(i)->getName()+")/dz", p.second[2]  );//

        of.printField();
      }
  }
}

DumpGradient::~DumpGradient() {
}

}


}



