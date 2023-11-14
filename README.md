# μephem

[![build](https://github.com/cjhoward/uephem/actions/workflows/build.yml/badge.svg)](https://github.com/cjhoward/uephem/actions/workflows/build.yml)

## Overview

μephem is a tiny ephemeris utility capable of generating orbital state vectors for the major bodies in our solar system. The μephem source is a single C-file with a few hundred lines of code, making it easily ported to other languages or modified to run on embedded systems.

## Usage

	uephem <file> <item> <t0> [<t1> <resolution>]

* `file` must be a JPL Planetary and Lunar Ephemerides binary ephemeris file. For an overview of the available ephemerides, see <https://ssd.jpl.nasa.gov/planets/eph_export.html>. Little-endian ephemeris files can be downloaded from <ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux>. Big-endian ephemeris files can be downloaded from <ftp://ssd.jpl.nasa.gov/pub/eph/planets/SunOS>. μephem is capable of reading both little and big-endian files, regardless of the host endianness. `.bsp` ephemeris files are **not** supported.
* `item` is the ID number of an ephemeris item to evaluate. See [the table below](#items) for more details.
* `t0` and `t1` are the start and end Julian dates. e.g. `2451545.0` corresponds to the J2000.0 epoch.
* `resolution` is the number of time points to evaluate between `t0` and `t1`.

### Items

| ID | Name                          | Output                                | Units     | Center |
|---:|:------------------------------|:--------------------------------------|:----------|:------:|
|  0 | Mercury                       | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  1 | Venus                         | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  2 | Earth-Moon barycenter         | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  3 | Mars                          | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  4 | Jupiter                       | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  5 | Saturn                        | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  6 | Uranus                        | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  7 | Neptune                       | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  8 | Pluto                         | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
|  9 | Moon                          | x,y,z,vx,vy,vz                        | km,km/d   | Earth  |
| 10 | Sun                           | x,y,z,vx,vy,vz                        | km,km/d   | SSB    |
| 11 | Earth nutation                | d(psi),d(epsilon),vd(psi),vd(epislon) | rad,rad/d | Earth  |
| 12 | Lunar mantle libration        | phi,theta,psi,vphi,vtheta,vpsi        | rad,rad/d | Moon   |
| 13 | Lunar mantle angular velocity | omega_x,omega_y,omega_z               | rad/d     | Moon   |
| 14 | TT-TDB                        | t                                     | s         | Earth  |

> Note: Items 11-14 may not be present in some ephemeris files.

## Configuration & Building

CMake is required to configure and build the application.

Run the following commands to configure and build a release executable:

	cmake -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build

## Ephemeris File Format

JPL ephemeris binary files consist of a series fixed-size records. Record size can vary between ephemeris files. A file's record size can be inferred through its file size along with the time information given in its header.

The first record in a file contains the header, the structure of which is described below:

	typedef struct
	{
		// Comment text.
		char comments[3][84];
		
		// First set of constant name labels.
		char cnames1[400][6];
		
		// JD start, JD end, record duration.
		double time[3];
		
		// Number of constants in the second record.
		int32_t nconst;
		
		// Constant values AU and EMRAT.
		double constants[2];
		
		// Record offset, coeffs per component, and coeff sets per record, for items 0-11.
		int32_t table1[12][3];
		
		// DE version number.
		int32_t denum;
		
		// Continuation of table1 for item 12.
		int32_t table2[3];
		
		// Second set of constant name labels.
		char cnames2[nconst - 400][6];
		
		// Continuation of table2 for items 13 and 14.
		int32_t table3[2][3];
		
	} de_header;

> Note: Some ephemeris files store additional constant names in the header between `table2` and `table3` if the limit of 400 constant names has been reached.
> 
> The second record is simply an array of `double` values for the constants named in the header. Each of the remaining records contain two `double` values: the Julian date start and end times of the record, followed by an array of `double` Chebyshev polynomial coefficients for each item in the ephemeris.

## License

μephem is licensed under the MIT license. See [LICENSE.md](./LICENSE.md) for details.
