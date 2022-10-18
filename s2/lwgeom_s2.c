#include "postgres.h"
#include "fmgr.h"


Datum s2_noop(PG_FUNCTION_ARGS);


/** Test library linking */
PG_FUNCTION_INFO_V1(s2_noop);
Datum s2_noop(PG_FUNCTION_ARGS)
{
	GSERIALIZED *in = PG_GETARG_GSERIALIZED_P(0);
	LWGEOM *lwgeom = lwgeom_from_gserialized(in);
	GSERIALIZED *out = geometry_serialize(lwgeom);
	PG_RETURN_POINTER(out);
}
