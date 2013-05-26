#pragma hdrstop
#include "rrRoadRunner.h"
#include "rrException.h"
#include "rrMinimizationData.h"
#include "rrc_minimization_api.h"
#include "rrc_utilities.h"   		//Support functions, not exposed as api functions and or data
//---------------------------------------------------------------------------

namespace rrc
{
using namespace std;
using namespace rr;

bool rrcCallConv addDoubleParameter(RRMinimizationDataHandle handle, const char* name, double value)
{
	try
    {
    	MinimizationData* data = castToMinimizationData(handle);
        data->addParameter(name, value);
        return true;
    }
    catch_bool_macro
}

bool rrcCallConv setMinimizationSelectionList(RRMinimizationDataHandle handle, RRStringArrayHandle listHandle)
{
	try
    {
    	MinimizationData* data = castToMinimizationData(handle);
//        StringList aList =

        return true;
    }
    catch_bool_macro
}


}