/*
 * EpipModelGenerator.cpp
 *
 * Created on: Jan 13, 2015 from GPUSimModelGenerator.cpp
 *
 * Author: JKM,
 */
#pragma hdrstop
#include "EpipModelGenerator.h"
#include "EpipExecutableModel.h"
// #include "ModelGeneratorContext.h"
#include "rrUtils.h"
#include "EpipException.h"
#include <rrLogger.h>
#include <Poco/Mutex.h>
#include <algorithm>

using rr::Logger;
using rr::getLogger;
using rr::ExecutableModel;
using rr::ModelGenerator;
using rr::Compiler;

namespace rr
{

namespace rrepip
{


EpipModelGenerator::EpipModelGenerator(const std::string &str)
: compilerStr(str)
{
    std::transform(compilerStr.begin(), compilerStr.end(), compilerStr.begin(), toupper);
    Log(Logger::LOG_TRACE) << __FUNC__;
}

EpipModelGenerator::~EpipModelGenerator()
{
    Log(Logger::LOG_TRACE) << __FUNC__;

}

bool EpipModelGenerator::setTemporaryDirectory(const std::string& path)
{
    return true;
}

std::string EpipModelGenerator::getTemporaryDirectory()
{
    return "not used";
}

ExecutableModel* EpipModelGenerator::createModel(const std::string& sbml,
        uint options)
{
    return new EpipExecutableModel(sbml, options);
}

Compiler* EpipModelGenerator::getCompiler()
{
    return NULL;
}

bool EpipModelGenerator::setCompiler(const std::string& compiler)
{
    return true;
}

} // namespace rrepip

} // namespace rr
