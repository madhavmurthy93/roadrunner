#pragma hdrstop
#include "rrRoadRunnerThread.h"

namespace rr
{
RoadRunnerThread::RoadRunnerThread(RoadRunner* rr) :
mRR(rr),
mIsDone(false)
{}


void RoadRunnerThread::setThreadName(const string& name)
{
	mThreadName = name;
}

void RoadRunnerThread::start()
{
	mThread.start(*this);
}

void RoadRunnerThread::run()	//This is called from runnable::start(this)
{
	worker();
}

void RoadRunnerThread::join()
{
	mThread.join();
}

void RoadRunnerThread::setRRInstance(RoadRunner* rr)
{
	mRR = rr;
}

RoadRunner* RoadRunnerThread::getRRInstance()
{
	return mRR;
}

bool RoadRunnerThread::isRunning()
{
	return mThread.isRunning();
}

}