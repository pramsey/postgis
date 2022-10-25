#include "postgres.h"
#include "fmgr.h"

#include "liblwgeom.h"
#include "lwgeom_pg.h"

Datum s2_noop(PG_FUNCTION_ARGS);


/*
 * This is required for builds against pgsql
 */
PG_MODULE_MAGIC;

/** Test library linking */
PG_FUNCTION_INFO_V1(s2_noop);
Datum s2_noop(PG_FUNCTION_ARGS)
{
	GSERIALIZED *in = PG_GETARG_GSERIALIZED_P(0);
	LWGEOM *lwgeom = lwgeom_from_gserialized(in);
	GSERIALIZED *out = geometry_serialize(lwgeom);
	PG_RETURN_POINTER(out);
}
