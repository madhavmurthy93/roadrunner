/*
 * EpipModelGenerator.h
 *
 * Created on: Jun 3, 2013
 *
 * Author: JKM
 */
#ifndef rrEpipModelGeneratorH
#define rrEpipModelGeneratorH

#include "ModelGenerator.h"

#if (__cplusplus >= 201103L) || defined(_MSC_VER)
#include <memory>
#include <unordered_map>
#define cxx11_ns std
#else
#include <tr1/memory>
#include <tr1/unordered_map>
#define cxx11_ns std::tr1
#endif

namespace rr
{

namespace rrepip
{

/**
 * General concepts:
 *
 * SBML defines chemical species.
 *
 * Floating Species: species whose value varies over time, i.e. they can 'float'
 * Boundary Species: boundary condition species, values are fixed to initial conditions.
 */
class RR_DECLSPEC EpipModelGenerator: public rr::ModelGenerator
{
public:
    EpipModelGenerator(const std::string &compilerStr);
    virtual ~EpipModelGenerator();

    /**
     * certain model generators, such as the compiler based ones
     * generate files such as shared libraries. This specifies the
     * location where they are stored.
     */
    virtual bool setTemporaryDirectory(const std::string& path);

    /**
     * certain model generators, such as the compiler based ones
     * generate files such as shared libraries. This specifies the
     * location where they are stored.
     */
    virtual std::string getTemporaryDirectory();

    /**
     * Create an executable model from an sbml string
     */
    virtual  rr::ExecutableModel *createModel(const std::string& sbml, uint options);


    /**
     * Get the compiler object that the model generator is using to
     * 'compile' sbml. Certain model generators may be interpreters, in this
     * case, the Compiler interface should still be sufficiently general to
     * manipulate interpreters as well.
     */
    virtual rr::Compiler *getCompiler();

    /**
     * No effect.
     */
    virtual bool setCompiler(const std::string& compiler);

private:
    std::string compilerStr;
};

} // namespace rrepip

} // namespace rr
#endif /* rrEpipModelGeneratorH */
