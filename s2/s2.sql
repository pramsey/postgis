
CREATE OR REPLACE FUNCTION s2_noop(geometry)
	RETURNS geometry
	AS '/Users/pramsey/Code/postgis-git/s2/postgis_s2-3.so', 's2_noop'
	LANGUAGE 'c' 
	IMMUTABLE STRICT 
	PARALLEL SAFE;

