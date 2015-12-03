
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

import java.io.*;
import java.util.*;

import mapper.GlobalOptions;
import mapper.Molecule;
import mapper.Reaction;

public class Mapper2000
{
	// Starting point of the program
	public static void main(String[] args) throws IOException
	{
		String arg = "";
		for (String s : args)
			arg += s + " ";

		if (arg.indexOf("-h") >= 0)
		{
			// print help
			System.out.println("Usage: java Mapper2000 [-algorithm -opts] [reactions]\n" +
					"Arguments:\n" +
					" search algorithm:\n" +
					"  -astar        - for astar search [default]\n" +
					"  -greedy       - for greedy search\n" +
					"  -dfs          - for dfs search\n" +
					"  -bpm          - for hungarian method\n" +
					"  -mcs          - for maximum common subgraph method\n" +
					"  -astar++ etc. - for astar using atom features\n" +
					" search params:\n" + 
					"  -one          - fetch only one solution\n" +
					"  -start n      - for starting at point 'n' [default=" + GlobalOptions.start + "]\n" +
					"  -end n        - for ending at point 'n-1'\n" +
					"  -bc           - count bond-changes as cost too\n" +
					"  -maxheap n    - max heap size, [default=" + GlobalOptions.MAX_HEAP + "]\n" +
					"  -k n          - greedy frequency [default=" + GlobalOptions.greedyFreq + "]\n" +
					"  -lb n         - accept solution with cost n\n" +
					"  -maxtime n    - max runtime in seconds [default=" + GlobalOptions.maxtime + "]\n" +
					"  -progress     - save progress of the lb/ub\n" +
					" dirs:\n" +
					"  -output       - output dir     [default=" + GlobalOptions.outputdir + "]\n" +
					"  -moldir       - mol dir        [default=" + GlobalOptions.molpath + "]\n" +
					"  -reacfile     - reaction list  [default=" + GlobalOptions.reacfile + "]\n" +
					"  -featdir      - feat directory [default=" + GlobalOptions.featpath + "]\n" +
					"  -overwrite    - overwrite results\n" +
					"  -nooutput     - no output\n" +
					"  -maxresults n - max number of results [default=" + GlobalOptions.maxresults + "]\n" +
					"  -rgraphs      - generate reaction graphs\n" + 
					" rest of the args - reaction id's\n" + 
					" without any reaction id's, do the automatical thing on whole kegg\n" +
					" if rest of the args contains a file, it's assumed to be the 'kegg_reactions.txt'\n\n" +
					"Copyright 2011 Markus Heinonen");
		
			System.exit(0);
		}		
		
		if (arg.indexOf("-dfs") >= 0)    { GlobalOptions.alg = "dfs"; GlobalOptions.plusplus = false; }
		if (arg.indexOf("-bpm") >= 0)    { GlobalOptions.alg = "bpm"; GlobalOptions.plusplus = false; }
		if (arg.indexOf("-greedy") >= 0) { GlobalOptions.alg = "greedy"; GlobalOptions.plusplus = false; }
		if (arg.indexOf("-mcs") >= 0)    { GlobalOptions.alg = "mcs"; GlobalOptions.plusplus = false; }
		if (arg.indexOf("-astar") >= 0)  { GlobalOptions.alg = "astar"; GlobalOptions.plusplus = false; }

		if (arg.indexOf("-one") >= 0) GlobalOptions.one = true;
		if (arg.indexOf("-overwrite") >= 0) GlobalOptions.overwrite = true;
		if (arg.indexOf("-nooutput") >= 0) GlobalOptions.noOutput = true;
		if (arg.indexOf("-bc") >= 0) GlobalOptions.bc = true;
		if (arg.indexOf("-progress") >= 0) GlobalOptions.progressinfo = true;
		if (arg.indexOf("-rgraphs") >= 0) GlobalOptions.rgraphs = true;
		if (arg.indexOf("-astar++") >= 0 || arg.indexOf("-greedy++") >= 0 || arg.indexOf("-bpm++") >= 0 || arg.indexOf("-dfs++") >= 0 || arg.indexOf("-mcs++") >= 0)
			GlobalOptions.plusplus = true;

		if (arg.indexOf("-start") >= 0)
			GlobalOptions.start = Mapper2000.getParam(arg, "-start");
		if (arg.indexOf("-end") >= 9)
			GlobalOptions.end = Mapper2000.getParam(arg, "-end");
		if (arg.indexOf("-lb") >= 0)
			GlobalOptions.lb = Mapper2000.getParam(arg, "-lb");
		if (arg.indexOf("-maxheap") >= 0)
			GlobalOptions.MAX_HEAP = Mapper2000.getLongParam(arg, "-maxheap");
		if (arg.indexOf("-maxresults") >= 0)
			GlobalOptions.maxresults = Mapper2000.getParam(arg, "-maxresults");
		if (arg.indexOf("-maxtime") >= 0)
			GlobalOptions.maxtime = Mapper2000.getParam(arg, "-maxtime");
		if (arg.indexOf("-k") >= 0)
			GlobalOptions.greedyFreq = Mapper2000.getParam(arg, "-k");
		if (arg.indexOf("-output") >= 0)
			GlobalOptions.outputdir = Mapper2000.getStrParam(arg, "-output");
		if (arg.indexOf("-moldir") >= 0)
			GlobalOptions.molpath = Mapper2000.getStrParam(arg, "-moldir");
		if (arg.indexOf("-reacfile") >= 0)
			GlobalOptions.reacfile = Mapper2000.getStrParam(arg, "-reacfile");
		if (arg.indexOf("-featdir") >= 0)
			GlobalOptions.featpath = Mapper2000.getStrParam(arg, "-featdir");
		

		// Get all reaction id's at the end of the arguments
		List<String> reac_ids = new ArrayList<String>();
		for (String s : args)
			if (s.charAt(0) == 'R')
				reac_ids.add(s);

		// no reaction names given, use batch mode for whole kegg
		if (reac_ids.size() == 0)
			GlobalOptions.batch = true;
		
		// Normal batch mode, do the thing against Kegg		
		if (GlobalOptions.batch)
		{
			try
			{
				runKeggReactions(GlobalOptions.reacfile, GlobalOptions.alg, GlobalOptions.start, GlobalOptions.end);
			} catch (Exception e)
			{
				System.out.println(e.getMessage());
				System.exit(1);
			}
		}
		else if (reac_ids.size() > 0) // otherwise, run against a give set of reactions
		{
			for (String s : reac_ids)
			{
				Reaction r = new Reaction(s);
				try
				{
					r.read();
				} catch (Exception e)
				{
					System.out.println(e.getMessage());
					continue;
				}

				System.gc();
				
				runReaction(r, GlobalOptions.alg);
			}
		}
	}

	/**
	 * Runs all kegg reactions, checks them from 'kegg_reactions.txt'
	 * 
	 */
	public static void runKeggReactions(String reac_file, String alg, int start, int end)
			throws Exception
	{
		Scanner sc = null;
		try
		{
			sc = new Scanner(new File(reac_file));
		} catch (IOException e)
		{
			throw new Exception("No kegg reaction file found: " + reac_file);
		}

		List<String> reac_strs = new ArrayList<String>();

		// throw away the first 'start-1' elements, default value for start=1, 
		// and thus we start from first element
		int i = 1;
		while (sc.hasNextLine() && i < start)
		{
			sc.nextLine();
			i++;
		}

		// read all remaining lines
		while (sc.hasNextLine() && i < end)
		{
			reac_strs.add(sc.nextLine());
			i++;
		}

		int j = 0;
		for (String line : reac_strs)
		{
			System.gc();
			
			System.out.print(j++ + start + ": ");

			String reac_id = line.split(" ")[0];

			Reaction r = new Reaction(reac_id);
			try
			{
				r.read();
			} catch (Exception e)
			{
				System.out.println(e.getMessage());
				continue;
			}

			runReaction(r, alg);
		}
	}

	public static void runReaction(Reaction r, String alg)
	{
		String output = "";
		
		if (GlobalOptions.outputdir.equals("./"))
		{
			output = "./" + alg;
			
			if (GlobalOptions.plusplus)
				output += "++";
			if (GlobalOptions.bc)
				output += "-bc";
			
			output += "-mappings/";
		}
		else
		{
			output = GlobalOptions.outputdir;
		}
		
		if (GlobalOptions.overwrite == false)
		{
			// Does the reaction file already exist
			File f;
			f = new File(output + r.getId() + ".txt");

			if (f.exists())
			{
				System.out.println("Reaction " + r.getId() + " already computed! Returning..");
				return;
			}
		}

		// Print something out
		System.out.println("Starting reaction " + r.getId());
		for (Molecule m : r.getSubs())
			System.out.println(" Substrate: " + m.getId() + ", " + m);
		for (Molecule m : r.getProds())
			System.out.println(" Product: " + m.getId() + ", " + m);

		long time = 0;
		
		try {
			time = r.computeMappings(alg);
		} catch (Exception e)
		{
			System.out.println(e.getMessage());
		}
		
		if (r.getMappings().size() == 0)
		{
//			System.out.println("Not enough memory");
			return;
		}

		int resultcount = r.getMappings().size();
		
		// there might be isomorphic mappings
		// -> test all mappings and only return truely different mappings (non-isomorphic)
		if (GlobalOptions.isomorphic)
		{
			r.contractIsomorphic();
			
			System.out.println(resultcount + " maps (" + r.getMappings().size() + " non-isomorphic), minimum cost " + r.getFirstMapping().costFunction() + ". Computation of maps took " + (time / 1000.0) + " secs.");
		}
		else
			System.out.println(resultcount + " maps, minimum cost " + r.getFirstMapping().costFunction() + ". Computation of maps took " + (time / 1000.0) + " secs.");
			
		
		if (!GlobalOptions.noOutput)
		{
			if (r.getMappings().size() != 0)
			{
				r.printMappingToFile(output, time);
				if (GlobalOptions.rgraphs)
					r.writeReactionGraphs(output);
			}
			else
				System.out.println("No mappings for reaction " + r.getId() + ". Sorry.");
		}
	}

	public static int getParam(String arg, String param)
	{
		int st = 0;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = Integer.parseInt(arg.substring(begin, stop).trim());
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}
	
	public static long getLongParam(String arg, String param)
	{
		long st = 0;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = Long.parseLong(arg.substring(begin, stop).trim());
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}

	public static String getStrParam(String arg, String param)
	{
		String st = null;
		try
		{
			int begin = arg.indexOf(param) + param.length() + 1;
			int stop = arg.indexOf(" ", begin);
			st = arg.substring(begin, stop).trim();
		} catch (NumberFormatException e)
		{
			System.out.println("Awkward argument for param " + arg + ", " + e);
			System.exit(0);
		}
		return st;
	}
}


