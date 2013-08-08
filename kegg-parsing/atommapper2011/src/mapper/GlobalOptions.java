//
//    Copyright 2011 Markus Heinonen 
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

package mapper;

public class GlobalOptions
{
	public static int start             = 1;
	public static int end               = Integer.MAX_VALUE;
	public static String alg            = "astar";
	public static boolean bc            = false;
	public static boolean overwrite     = false;
	public static boolean isomorphic    = true;
	public static boolean noOutput      = false;
	public static boolean batch         = false;
	public static boolean rgraphs       = false;
	public static long MAX_HEAP         = Integer.MAX_VALUE;
	public static int maxresults        = Integer.MAX_VALUE;
	public static int greedyFreq        = 1;
	public static boolean plusplus      = true;
	public static String outputdir      = "./";
	public static int maxtime           = Integer.MAX_VALUE;
	public static boolean progressinfo  = false;
	public static int lb                = 0;
	public static boolean one           = false;
	public static String molpath        = "/group/home/icomic/data/kegg/ligand/2010.07.01/mol/";
	public static String reacfile       = "../kegg-reactions.txt";
	public static String featpath       = "/home/group/icomic/data/kegg/mol-features-2010.07.01/";
}

