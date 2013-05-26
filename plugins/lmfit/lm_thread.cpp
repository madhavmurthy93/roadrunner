#pragma hdrstop
#include "rrRoadRunnerData.h"
#include "rrLogger.h"
#include "lm_thread.h"
#include "rrNoise.h"
#include "lm.h"
#include "lmfit-3.5/lib/lmcurve.h"
#include "../../Wrappers/C/rrc_api.h"
#include "../../Wrappers/C/rrc_support.h"

//Global roadrunner... because lm fit don't allow void* to be passed in...
RoadRunner* rri;
Parameters  gParas;

//---------------------------------------------------------------------------
using namespace rr;
LMFitThread::LMFitThread(LM& host)
:
threadEnterCB(NULL),
threadExitCB(NULL),
mUserData(NULL),
mTheHost(host),
mMinData(mTheHost.getMinimizationData())
{}

bool LMFitThread::isRuning()
{
	return mThread.isRunning();
}

void LMFitThread::assignCallBacks(ThreadCB fn1, ThreadCB fn2, void* userData)
{
	threadEnterCB 	= fn1;
	threadExitCB  	= fn2;
    mUserData 		= userData;
}

void LMFitThread::start()
{
    if(mThread.isRunning())
    {
    	Log(lError)<<"Tried to start an already working thread!";
        return;
    }

	mThread.start(*this);	//Will spawn the run function below in a thread
}

double f(double t, const double* paras);

void LMFitThread::run()
{
	if(threadEnterCB)
    {
		threadEnterCB(mUserData);	//Tell anyone who wants to know
    }

   	if(rri)
    {
    	delete rri;
    }

    rri = new RoadRunner;
    rri->loadSBML(mTheHost.mSBML.getValue(), false);

    Log(lInfo)<<"The following parameters are to be minimized";
    gParas =  mMinData.getParameters();

    for(int i = 0; i < gParas.count(); i++)
    {
    	Log(lInfo)<<gParas[i]->getName()<<" with initial value: "<<gParas[i]->getValueAsString();
    }

    RoadRunnerData inputData = mMinData.getInputData();
    Log(lInfo)<<"The following data is to be used for fitting";
    Log(lInfo)<<inputData;

	//Create data vector appropriate for lmcurve_fit
    int 	n_par = gParas.count();
    double *par   = new double[n_par];

    for(int i = 0; i < n_par; i++)
    {
    	string val 	= gParas[i]->getValueAsString();
    	par[i] 		= toDouble(val);
    }

    int count = inputData.rSize();
    double *time 	= new double[count];
    double *y 		= new double[count];

    //populate time and y
    for(int i = 0; i < count; i++)
    {
    	time[i] = inputData(i,0);
    	y[i] 	= inputData(i,1);
    }

    //Aux parameters
    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    control.printflags = 3;

    rri->setTimeCourseSelectionList("S1");

    lmcurve_fit(    n_par,
    			    par,
                    count,
                    time,
                    y,
                    f,
                    &control,
                    &status
    			);

    Log(lInfo)<<"The LM algorithm finished with the following status: "<<lm_infmsg[status.info];
    Log(lInfo)<<"After "<<status.nfev<<" evaluations.";
    for(int i = 0; i < gParas.count(); i++)
    {
	    Log(lInfo)<<"Parameter "<< gParas[i]->getName() <<" is "<<par[i];
    }
    delete [] y;
    delete [] time;

    rri->reset();
    //Calculate model data
    RoadRunnerData modelData;
    RoadRunnerData residualsData;
    modelData.reSize(inputData.rSize(), inputData.cSize());
    residualsData.reSize(inputData.rSize(), inputData.cSize());
	for(int row = 0; row < inputData.rSize(); row++)
    {
        //Time....
		modelData(row, 0) = inputData(row,0);
		residualsData(row, 0) = inputData(row,0);

        //Y
        modelData(row, 1) = f(inputData(row,0), par);
        residualsData(row, 1) = inputData(row, 1) - modelData(row, 1);
    }

    //Populate minDataObject
	mMinData.setModelData(modelData);
	mMinData.setResidualsData(residualsData);

	if(threadExitCB)
    {
		threadExitCB(mUserData);
    }
}

double f(double time, const double* paras)
{
	double timeStep = 1.e-7;

	//set model parameters, contained in gParas
    for(int i = 0; i < gParas.count(); i++)
    {
    	Parameter<double>* aPar = (Parameter<double>*) gParas[i];
    	setValue(rri, aPar->getName().c_str(), paras[i]);
		aPar->setValue(paras[i]);
    }

    if(time != 0)
    {
		//propagate the model to time t
    	rri->simulateEx(0, time, 1);
    }

    //Get the values at current time...
	vector<SelectionRecord> sel = rri->getSelectionList();

    //for each specie, sum up current values
    double value = 0;
    for(int i = 0; i < sel.size(); i++)
    {
    	//value += rri->getValueForRecord(sel[i]);
		value = rri->getValueForRecord(sel[i]);
    }

    return value;
}

