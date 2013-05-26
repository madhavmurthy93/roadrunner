#ifndef rrMinimizationDataH
#define rrMinimizationDataH
#include <ostream>
#include "rrObject.h"
#include "rrStringlist.h"
#include "rrRoadRunnerData.h"
#include "rrParameter.h"
#include "rrParameters.h"
//---------------------------------------------------------------------------

namespace rr
{
using std::ostream;


class RR_DECLSPEC MinimizationData : public rrObject
{
	protected:
		RoadRunnerData                 	mExperimentalData;			//Observed data
		RoadRunnerData                 	mModelData;					//Observed data
		RoadRunnerData                 	mResidualsData;				//Observed data
		stringstream			       	mReport;
        Parameters				    	mParameters;				//Parameters to fit
        StringList						mSelectionList;


    public:
					                   	MinimizationData();
					                   ~MinimizationData();
					                   	MinimizationData(const MinimizationData& data);
		MinimizationData&				operator=(MinimizationData& rhs);
		void					       	init();

        void							addParameter(const string& name, const double& value);
        void							addParameter(const string& name, const int& value);
        void							setSelectionList(const StringList& selList);
        void							setInputData(RoadRunnerData& data);
        void							setModelData(RoadRunnerData& data);
        void							setResidualsData(RoadRunnerData& data);

        RoadRunnerData					getInputData();
        RoadRunnerData					getModelData();
        RoadRunnerData					getResidualsData();

        RoadRunnerData&					getInputDataReference();
        RoadRunnerData&					getModelDataReference();
        RoadRunnerData&					getResidualsDataReference();

        ostream&	   			       	operator<<(const string& str);
        string					       	getReport();
		Parameters					  	getParameters();
        bool							reset();

};

template<>
std::string Parameter< MinimizationData >::getType() const
{
    return "Pointer To Minimization Result";
}

template<>
string Parameter< MinimizationData >::getValueAsString() const
{
    throw("You can't get the value of this structure as string.. :(");
}

template<>
void Parameter< MinimizationData >::setValueFromString(const string& val)
{
	//We can't setup this data structure from a string... :(
    throw("You can't set the value of this structure from a string.. :(");
}

template<>
void Parameter< MinimizationData >::setValue(MinimizationData* val)
{
	mValue = *(val);
}

}
#endif