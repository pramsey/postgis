
CREATE OR REPLACE FUNCTION s2_noop(geometry)
	RETURNS geometry
	AS 'MODULE_PATHNAME', 's2_noop'
	LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
	_COST_DEFAULT;
