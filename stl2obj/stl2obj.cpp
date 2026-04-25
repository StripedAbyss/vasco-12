// stl2obj converts an STL CAD file to OBJ format.

// Copyright (c) 2017 Amir Baserinia

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdio>
#include <cstdlib>

#include "vectornd.h"
#include "geometry.h"
#include "importstl.h"
#include "exportobj.h"

// The name of this program
static const char* PROGRAM_NAME = "stl2obj";

// author
static const char* AUTHOR = "Amir Baserinia";

// usage help
void usage (int status)
{
    printf ("Usage: %s [OPTION]... [FILE]...\n", PROGRAM_NAME);
    printf ("Converts CAD STL models to OBJ format.\n");
    printf (
        "Options:\n"
        "  -m, --merge-vertices     merge vertices\n"
        "  -f, --fill-holes         file holes in surface\n"
        "  -s, --stich-cureves      stick curves between surfaces\n"
        "  -t, --tolerance          merge tolerance\n");
    printf (
        "Examples:\n"
        "  %s input.stl output.obj  convert input from STL to OBJ and write "
        "to output.\n", PROGRAM_NAME);
    exit (status);
}

// version information
void version ()
{
    printf ("%s converts an STL CAD file to OBJ format.\n", PROGRAM_NAME);
    printf ("Copyright (c) 2017 %s\n", AUTHOR);
}

int main (int argc, char **argv)
{

//  create a geometry tesselation object
    Geometry tessel;

//  fill up the tesselation object with STL data (load STL)
    tessel.visit (ImportSTL ("C:\\Users\\Administrator\\Downloads\\stl2obj-master\\example\\Fidgit.stl"));

//  write down the tesselation object into OBJ file (save OBJ)
    tessel.visit (ExportOBJ ("C:\\Users\\Administrator\\Downloads\\stl2obj-master\\example\\Fidgit_2.obj"));

}

