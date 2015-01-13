/*
 * ModelGenerator.cpp
 *
 *  Created on: May 20, 2013
 *      Author: andy
 */
#pragma hdrstop
#include "ModelGenerator.h"


#if defined(BUILD_EPIPHANY)
#include "epiphany/EpipModelGenerator.h"
#endif

#if defined(BUILD_LLVM)
#include "llvm/LLVMModelGenerator.h"
#endif

#if defined(BUILD_LEGACY_C)
#include "c/rrCModelGenerator.h"
#endif

#include "rrLogger.h"
#include <string>
#include <algorithm>

using namespace std;
namespace rr {


std::string Compiler::getDefaultCompiler()
{
#if defined(BUILD_LLVM)
    return "LLVM";
#else
    #if defined(_WIN32)
        return joinPath("..", "compilers", "tcc", "tcc.exe");
    #else
        // the default compiler on Unix systems is 'cc', the standard enviornment
        // for the default compiler is 'CC'.
        return getenv("CC") ? getenv("CC") : "gcc";
    #endif
#endif
}

ModelGenerator* ModelGenerator::New(const string& compiler, const string& tempFolder,
            const string& supportCodeFolder)
{
    Log(Logger::LOG_INFORMATION) << "createing model generator, compiler: \"" << compiler << "\"";

    string ucomp = compiler;
    std::transform(ucomp.begin(), ucomp.end(),ucomp.begin(), ::toupper);

#if defined(BUILD_EPIPHANY)
    if (ucomp == "EPIPHANY")
    {
        Log(Logger::LOG_INFORMATION) << "Creating Epiphany based model generator.";
        return new rrepip::EpipModelGenerator(compiler);
    }
# endif

#if defined(BUILD_LLVM) && !defined(BUILD_LEGACY_C)

    Log(Logger::LOG_INFORMATION) << "Creating LLVM based model generator.";
    return new rrllvm::LLVMModelGenerator(compiler);

#endif

#if defined(BUILD_LLVM) && defined(BUILD_LEGACY_C)

    if (ucomp == "LLVM")
    {
        Log(Logger::LOG_INFORMATION) << "Creating LLVM based model generator.";
        return new rrllvm::LLVMModelGenerator(compiler);
    }

    Log(Logger::LOG_INFORMATION) << "Creating C based model generator using " << compiler << " compiler.";

    return new CModelGenerator(tempFolder, supportCodeFolder, compiler);
#endif

#if !defined(BUILD_LLVM) && defined(BUILD_LEGACY_C)

    Log(Logger::LOG_INFORMATION) << "Creating C based model generator using " << compiler << " compiler.";

    // default (for now...), the old C code generating model generator.
    return new CModelGenerator(tempFolder, supportCodeFolder, compiler);
#endif

#if !defined(BUILD_LLVM) && !defined(BUILD_LEGACY_C)
#error Must built at least one ModelGenerator backend, either BUILD_LLVM or BUILD_LEGACY_C
#endif


}

} /* namespace rr */
