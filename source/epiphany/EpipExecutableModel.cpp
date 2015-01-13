/*
 * EpipExecutableModel.cpp
 *
 * Author: JKM
 */
#pragma hdrstop
#include "EpipExecutableModel.h"
#include "rrSparse.h"
#include "rrLogger.h"
#include "rrException.h"
#include "rrStringUtils.h"
#include "rrConfig.h"
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <iostream>

/*
#if RR_GPUSim_USE_LLVM_MODEL
#include "llvm/LLVMModelGenerator.h"
#endif
*/

using rr::Logger;
using rr::getLogger;
using rr::LoggingBuffer;
using rr::SelectionRecord;
using rr::EventListener;
using rr::EventListenerPtr;
using rr::EventListenerException;
using rr::Config;

#if defined (_WIN32)
#define isnan _isnan
#else
#define isnan std::isnan
#endif

template <typename numeric_type>
static void dump_array(std::ostream &os, int n, const numeric_type *p)
{
    if (p)
    {
        os << setiosflags(std::ios::floatfield) << std::setprecision(8);
        os << '[';
        for (int i = 0; i < n; ++i)
        {
            os << p[i];
            if (i < n - 1)
            {
                os << ", ";
            }
        }
        os << ']' << std::endl;
    }
    else
    {
        os << "NULL" << std::endl;
    }
}

namespace rr
{

namespace rrepip
{

// -- GPUEntryPoint --
/*
GPUEntryPoint::GPUEntryPoint(void* sym, Precision p) {
    switch (p) {
        case Precision::Single:
            symsp_ = (EntryPointSigSP)sym;
            break;
        case Precision::Double:
            symdp_ = (EntryPointSigDP)sym;
            break;
//         default: // Shut up clang
//             assert(0 && "Should not happen");
    }
}

void GPUEntryPoint::operator()(int m, float* t, float* v) {
    if (symsp_)
        symsp_(m, t, v);
    else
        throw_Epip_exception("Wrong signature for this module");
}

void GPUEntryPoint::operator()(int m, double* t, double* v) {
    if (symdp_)
        symdp_(m, t, v);
    else
        throw_Epip_exception("Wrong signature for this module");
}
*/
  /**
 * checks if the bitfield value has all the flags
 * in type (equiv to checkExact but with a more accurate
 * name)
 */

inline bool checkBitfieldSubset(uint32_t type, uint32_t value) {
    return (value & type) == type;
}

EpipExecutableModel::EpipExecutableModel(std::string const &sbml, unsigned loadSBMLOptions) {
    // generate a hash of the SBML (TODO: find a faster way to get a UUID for the model)
    /*
    sbmlhash_.reset(new Hashval(Hash::me(sbml)));
    generator_->setPrecision(dom::CudaGenerator::Precision::Double);
#if RR_Epip_USE_LLVM_MODEL
    rrllvm::LLVMModelGenerator modelGenerator{Compiler::getDefaultCompiler()};
    llvmmodel_.reset(dynamic_cast<rrllvm::LLVMExecutableModel*>(modelGenerator.createModel(sbml, loadSBMLOptions)));
    assert(llvmmodel_ && "Internal error (no model)");
#endif
    */
    std::cout << "Hello World!!" << std::endl; 
}

EpipExecutableModel::~EpipExecutableModel() {
    Log(Logger::LOG_DEBUG) << __FUNC__;
}

/*
void EpipExecutableModel::generateModel() {
    generator_->generate(*this);
}

GPUEntryPoint EpipExecutableModel::getEntryPoint() {
    return generator_->getEntryPoint();
}
*/

string EpipExecutableModel::getModelName() {
    /*
    if(!getModel())
        throw_Epip_exception("No model loaded");
    return getModel()->getName();
    */
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::setTime(double time)
{
    throw_Epip_exception("not supported");
}

double EpipExecutableModel::getTime()
{
    throw_Epip_exception("not supported");
}


int EpipExecutableModel::getNumIndFloatingSpecies()
{
    throw_Epip_exception("not supported");
}


int EpipExecutableModel::getNumDepFloatingSpecies()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumFloatingSpecies()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumBoundarySpecies()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumGlobalParameters()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumCompartments()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumReactions()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumLocalParameters(int reactionId)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::convertToAmounts()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::computeConservedTotals()
{
    throw_Epip_exception("not supported");
}


int EpipExecutableModel::getFloatingSpeciesConcentrations(int len, int const *indx, double *values) {
    Log(Logger::LOG_TRACE) << "EpipExecutableModel::getFloatingSpeciesConcentrations: len = " << len << ", indx = " << indx << ", values = " << values;
/*
#if RR_Epip_USE_LLVM_MODEL
//     throw_Epip_exception("EpipExecutableModel::getFloatingSpeciesConcentrations");
    checkLLVMModel();
    return llvmmodel_->getFloatingSpeciesConcentrations(len,  indx,  values);
#else
    for (int k=0; k<len; ++k) {
        values[indx[k]] = 123.;
    }
    return len;
//     throw_Epip_exception("not supported");
#endif
*/
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::getRateRuleValues(double *rateRuleValues)
{
    throw_Epip_exception("not supported");
}


void EpipExecutableModel::convertToConcentrations()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::updateDependentSpeciesValues()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::computeAllRatesOfChange()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::getStateVectorRate(double time, const double *y, double *dydt) {
    if (!y && !dydt)
        // does nothing
        return;
    throw_Epip_exception("not supported");
}

double EpipExecutableModel::getFloatingSpeciesAmountRate(int index,
           const double *reactionRates)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::testConstraints()
{
}

std::string EpipExecutableModel::getInfo()
{
    std::stringstream stream;

    double *tmp;

    int nFloat = getNumFloatingSpecies();
    int nBound = getNumBoundarySpecies();
    int nComp = getNumCompartments();
    int nGlobalParam = getNumGlobalParameters();
    int nEvents = getNumEvents();
    int nReactions = getNumReactions();

    stream << "* Calculated Values *" << std::endl;

    tmp = new double[nFloat];
    getFloatingSpeciesAmounts(nFloat, 0, tmp);
    stream << "FloatingSpeciesAmounts:" << std::endl;
    dump_array(stream, nFloat, tmp);

    getFloatingSpeciesConcentrations(nFloat, 0, tmp);
    stream << "FloatingSpeciesConcentrations:" << std::endl;
    dump_array(stream, nFloat, tmp);

    this->getFloatingSpeciesInitConcentrations(nFloat, 0, tmp);
    stream << "FloatingSpeciesInitConcentrations:" << std::endl;
    dump_array(stream, nFloat, tmp);
    delete[] tmp;

    tmp = new double[nReactions];
    getReactionRates(nReactions, 0, tmp);
    stream << "Reaction Rates:" << std::endl;
    dump_array(stream, nReactions, tmp);
    delete[] tmp;

    tmp = new double[nBound];
    getBoundarySpeciesAmounts(nBound, 0, tmp);
    stream << "BoundarySpeciesAmounts:" << std::endl;
    dump_array(stream, nBound, tmp);

    getBoundarySpeciesConcentrations(nBound, 0, tmp);
    stream << "BoundarySpeciesConcentrations:" << std::endl;
    dump_array(stream, nBound, tmp);
    delete[] tmp;

    tmp = new double[nComp];
    getCompartmentVolumes(nComp, 0, tmp);
    stream << "CompartmentVolumes:" << std::endl;
    dump_array(stream, nComp, tmp);

    this->getCompartmentInitVolumes(nComp, 0, tmp);
    stream << "CompartmentInitVolumes:" << std::endl;
    dump_array(stream, nComp, tmp);
    delete[] tmp;

    tmp = new double[nGlobalParam];
    getGlobalParameterValues(nGlobalParam, 0, tmp);
    stream << "GlobalParameters:" << std::endl;
    dump_array(stream, nGlobalParam, tmp);
    delete[] tmp;

    tmp = new double[nGlobalParam];
    getGlobalParameterValues(nGlobalParam, 0, tmp);
    stream << "GlobalParameters:" << std::endl;
    dump_array(stream, nGlobalParam, tmp);
    delete[] tmp;

    unsigned char *tmpEvents = new unsigned char[nEvents];
    getEventTriggers(nEvents, 0, tmpEvents);
    stream << "Events Trigger Status:" << std::endl;
    dump_array(stream, nEvents, (bool*)tmpEvents);
    delete[] tmpEvents;

//     stream << *modelData;

    return stream.str();
}

int EpipExecutableModel::getFloatingSpeciesIndex(const string& id) {
    //return getFloatingSpeciesId(id)->getIndex();
    // temp return val = 1
    return 1;
}

string EpipExecutableModel::getFloatingSpeciesId(int index)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getBoundarySpeciesIndex(const string& id)
{
    throw_Epip_exception("not supported");
}

string EpipExecutableModel::getBoundarySpeciesId(int indx)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getGlobalParameterIndex(const string& id)
{
    throw_Epip_exception("not supported");
}

string EpipExecutableModel::getGlobalParameterId(int id)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getCompartmentIndex(const string& id)
{
    throw_Epip_exception("not supported");
}

string EpipExecutableModel::getCompartmentId(int id)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getReactionIndex(const string& id)
{
    throw_Epip_exception("not supported");
}

string EpipExecutableModel::getReactionId(int id)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::evalInitialConditions()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::reset()
{
    uint options = rr::Config::getInt(rr::Config::MODEL_RESET);
    Log(Logger::LOG_DEBUG) << "calling reset with default values: " << options;
    reset(options);
}

void EpipExecutableModel::reset(int options)
{
//     throw_Epip_exception("not supported");
    // resets the initial conditions
    if ((options & SelectionRecord::INITIAL)) {
        Log(Logger::LOG_INFORMATION) << "resetting init conditions";
//         evalInitialConditions();
    }

#if RR_Epip_USE_LLVM_MODEL
    checkLLVMModel();
    llvmmodel_->reset(options);
#endif
}

bool EpipExecutableModel::getConservedSumChanged()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::setConservedSumChanged(bool val)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getStateVector(double* stateVector)
{
#if RR_Epip_USE_LLVM_MODEL
    checkLLVMModel();
    return llvmmodel_->getStateVector(stateVector);
#else
    Log(Logger::LOG_TRACE) << "EpipExecutableModel::getStateVector: len of state vector = " << getNumIndFloatingSpecies();
    if (!stateVector)
        return getNumIndFloatingSpecies();
    throw_Epip_exception("getStateVector (non-null) not supported");
#endif
}

int EpipExecutableModel::setStateVector(const double* stateVector)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::print(std::ostream &stream)
{
    stream << "EpipExecutableModel" << std::endl;
    stream << getInfo();
}

void EpipExecutableModel::getIds(int types, std::list<std::string> &ids) {
    // EpipModel::getIds(types,ids);
}

int EpipExecutableModel::getSupportedIdTypes() {
    return SelectionRecord::TIME |
        SelectionRecord::BOUNDARY_CONCENTRATION |
        SelectionRecord::FLOATING_CONCENTRATION |
        SelectionRecord::REACTION_RATE |
        SelectionRecord::FLOATING_AMOUNT_RATE |
        SelectionRecord::FLOATING_CONCENTRATION_RATE |
        SelectionRecord::COMPARTMENT |
        SelectionRecord::GLOBAL_PARAMETER |
        SelectionRecord::FLOATING_AMOUNT |
        SelectionRecord::BOUNDARY_AMOUNT |
        SelectionRecord::INITIAL_AMOUNT |
        SelectionRecord::INITIAL_CONCENTRATION |
        SelectionRecord::STOICHIOMETRY;
}

double EpipExecutableModel::getValue(const std::string& id)
{
    throw_Epip_exception("not supported");
}


const rr::SelectionRecord& EpipExecutableModel::getSelection(const std::string& str)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::setValue(const std::string& id, double value)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getFloatingSpeciesConcentrationRates(int len,
        const int* indx, double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setBoundarySpeciesAmounts(int len, const int* indx,
        const double* values)
{
    throw_Epip_exception("not supported");
}

std::string EpipExecutableModel::getStateVectorId(int index)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::evalReactionRates()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumRateRules()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getFloatingSpeciesAmounts(int len, const int* indx,
        double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setFloatingSpeciesConcentrations(int len,
        const int* indx, const double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getBoundarySpeciesAmounts(int len, const int* indx,
        double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getBoundarySpeciesConcentrations(int len,
        const int* indx, double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setBoundarySpeciesConcentrations(int len,
        const int* indx, const double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getGlobalParameterValues(int len, const int* indx,
        double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setGlobalParameterValues(int len, const int* indx,
        const double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getCompartmentVolumes(int len, const int* indx,
        double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getReactionRates(int len, const int* indx,
        double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getNumConservedMoieties()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getConservedMoietyIndex(const string& name)
{
    throw_Epip_exception("not supported");
}

string EpipExecutableModel::getConservedMoietyId(int index)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getConservedMoietyValues(int len, const int* indx,
        double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setConservedMoietyValues(int len, const int* indx,
        const double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getFloatingSpeciesAmountRates(int len,
        int const *indx, double *values)
{
    throw_Epip_exception("not supported");
}


int EpipExecutableModel::setFloatingSpeciesAmounts(int len, int const *indx,
        const double *values)
{
    throw_Epip_exception("not supported");
}


int EpipExecutableModel::setCompartmentVolumes(int len, const int* indx,
        const double* values)
{
    throw_Epip_exception("not supported");
}



double EpipExecutableModel::getStoichiometry(int speciesIndex, int reactionIndex)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getStoichiometryMatrix(int* pRows, int* pCols,
        double** pData)
{
    throw_Epip_exception("not supported");
}



/******************************* Initial Conditions Section *******************/
#if (1) /**********************************************************************/
/******************************************************************************/


int EpipExecutableModel::setFloatingSpeciesInitConcentrations(int len,
        const int* indx, const double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getFloatingSpeciesInitConcentrations(int len,
        const int* indx, double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setFloatingSpeciesInitAmounts(int len, int const *indx,
            double const *values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getFloatingSpeciesInitAmounts(int len, int const *indx,
                double *values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setCompartmentInitVolumes(int len, const int *indx,
            double const *values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getCompartmentInitVolumes(int len, const int *indx,
                double *values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::setGlobalParameterInitValues(int len, const int* indx,
        const double* values)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getGlobalParameterInitValues(int len, const int *indx,
                double *values)
{
    throw_Epip_exception("not supported");
}

/******************************* End Initial Conditions Section ***************/
#endif /***********************************************************************/
/******************************************************************************/

/******************************* Events Section *******************************/
#if (1) /**********************************************************************/
/******************************************************************************/

int EpipExecutableModel::getNumEvents() {
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getEventTriggers(int len, const int *indx, unsigned char *values)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::applyEvents(double timeEnd, const unsigned char* previousEventStatus,
	    const double *initialState, double* finalState)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::applyPendingEvents(const double *stateVector, double timeEnd, double tout)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::getEventRoots(double time, const double* y, double* gdot)
{
    throw_Epip_exception("not supported");
}

double EpipExecutableModel::getNextPendingEventTime(bool pop)
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getPendingEventSize()
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::resetEvents()
{
    throw_Epip_exception("not supported");
}

int EpipExecutableModel::getEventIndex(const std::string& eid)
{
    throw_Epip_exception("not supported");
}

std::string EpipExecutableModel::getEventId(int index)
{
    throw_Epip_exception("not supported");
}

void EpipExecutableModel::setEventListener(int index, rr::EventListenerPtr eventHandler)
{
    throw_Epip_exception("not supported");
}

rr::EventListenerPtr EpipExecutableModel::getEventListener(int index)
{
    throw_Epip_exception("not supported");
}

/******************************* Events Section *******************************/
  #endif /**********************************************************************/
/******************************************************************************/

} // namespace rrepip

} // namespace rr
