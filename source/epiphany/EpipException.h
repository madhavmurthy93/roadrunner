/*
 * EpipException.h
 *
 *  Created on: Jun 30, 2013
 *      Author: JKM
 */

#ifndef RREpipEXCEPTION_H_
#define RREpipEXCEPTION_H_

#include <stdexcept>
#include "rrLogger.h"
#include "rrOSSpecifics.h"

namespace rr
{

namespace rrepip
{

class EpipException: public std::runtime_error
{
public:
    explicit EpipException(const std::string& what) :
            std::runtime_error(what)
    {
    }

    explicit EpipException(const std::string& what, const std::string &where) :
            std::runtime_error(what + ", at " + where)
    {
    }
};

#define throw_Epip_exception(what) \
        {  \
            Log(rr::Logger::LOG_INFORMATION) << "EpipException, what: " \
                << what << ", where: " << __FUNC__; \
                throw EpipException(what, __FUNC__); \
        }


} // namespace rrepip

} // namespace rr

#endif /* RREpipEXCEPTION_H_ */
