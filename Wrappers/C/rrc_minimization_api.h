#ifndef rrc_minimization_apiH
#define rrc_minimization_apiH
#include "rrc_exporter.h"
#include "rrc_types.h"
//---------------------------------------------------------------------------

#if defined(__cplusplus)
namespace rrc
{
extern "C"
{
#endif

/*!
 \brief Add parameter to fit to minimization data structure
 \param[in] RRMinimizationDataHandle
 \param[in] char* name
 \param[in] double value
 \return Returns true if sucessful, false otherwise
 \ingroup Minimization
*/
C_DECL_SPEC bool rrcCallConv addDoubleParameter(RRMinimizationDataHandle handle, const char* name, double value);

/*!
 \brief Set minimization objects selection list
 \param[in] RRMinimizationDataHandle
 \param[in] RRStringArrayHandle
 \return Returns true if sucessful, false otherwise
 \ingroup Minimization
*/
C_DECL_SPEC bool rrcCallConv setMinimizationSelectionList(RRMinimizationDataHandle handle, RRStringArrayHandle listHandle);


#if defined(__cplusplus)
}	//Extern "C"

}	//rrc namespace
#endif

#endif
