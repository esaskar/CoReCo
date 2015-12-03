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

import java.io.*;
import java.util.*;


public class Reaction
{
	// e.g. R04322
	private String id = "no id";
	private String reaction = "no reaction";
	private int molsNum;
	private List<Molecule> substrates = new ArrayList<Molecule>();
	private List<Molecule> products = new ArrayList<Molecule>();
	private List<Mapping> mappings = new ArrayList<Mapping>();
	
	private Collection<Atom> subatoms = new HashSet<Atom>();
	private Collection<Atom> prodatoms = new HashSet<Atom>();
	private Collection<Bond> subbonds = new HashSet<Bond>();
	private Collection<Bond> prodbonds = new HashSet<Bond>();
	
	private HashMap<Atom, HashMap<Atom, Integer>> distances = new HashMap<Atom, HashMap<Atom,Integer>>();
	private Collection<Bond> removedBonds = new HashSet<Bond>();
	private Collection<Atom> removedAtoms = new HashSet<Atom>();
	
	private double[][][] morgan;
	private double[][][] wiener;
	private double[][][] ad;
	private double[][][] bd;

	private double[][] featcost;
	private double[][] similarity;
	
	
	public Reaction(String id)
	{
//		Molecule.resetMolCounter();
//		Atom.resetReactionalAtomCounter();
		this.id = id;
	}
	
	public Reaction(List<Molecule> subs, List<Molecule> prods)
	{
		this.substrates = subs;
		this.products = prods;
		this.molsNum = subs.size() + prods.size();
	}

	public Reaction(List<Molecule> subs, List<Molecule> prods,
			List<Mapping> maps)
	{
		this(subs, prods);
		this.mappings = maps;
	}

	public Reaction(List<Molecule> subs, List<Molecule> prods, Mapping m)
	{
		this(subs, prods);
		this.mappings.add(m);
	}

	public Mapping getFirstMapping()
	{
		return this.mappings.get(0);
	}
	
	public List<Mapping> getMappings()
	{
		return this.mappings;
	}

	public String getId()
	{
		return this.id;
	}
	
	public boolean equals(Object ob)
	{
		return id == ((Reaction)ob).getId();
	}
	
	public int hashCode()
	{
		return this.id.hashCode();
	}

	// returns time spent on computation
	public long computeMappings(String alg) throws Exception
	{
		Timer t = new Timer();
		
		if (GlobalOptions.plusplus)
			this.constructFeatureTables();
		
		computeAtomStringIndex();
		computeBondStringIndex();
		
//		int ub = computeInitialUB();
//		int lb = computeInitialLB();
		
		if (alg == "greedy")
		{
			t.start();
			Mapping m = Greedy.greedy(this);
			this.mappings.add(m);
			t.stop();
		}
		else if (alg == "dfs")
		{
			t.start();
			this.mappings = new ArrayList<Mapping>( DFS.runDFS(this) );
			t.stop();
		}
		else if (alg == "astar")
		{
			t.start();
			
			Collection<Mapping> temp = null;
			
			try {
				temp = Astar.astarFixedOrder(this);
			}
			catch (Exception e)
			{
				System.out.println(e.getMessage());
			}
				
			if (temp == null)
				return 0;
			
			this.mappings = new ArrayList<Mapping>( temp );
			t.stop();
		}
		else if (alg == "bpm")
		{
			t.start();
			this.mappings = Hungarian.match(this);
			t.stop();
		}
		else if (alg == "mcs")
		{
			t.start();
			this.mappings = MCS.mcs(this);
			t.stop();
		}
		
		
		return t.getTime();
	}
	
	// computes the initial upper bound based on the reaction
	// we have to make at maximum
	public int computeInitialUB()
	{
//		return Integer.MAX_VALUE;
		return getSubsBonds().size() + getProdsBonds().size();
	}
	
	// initial lower bound
	// we have to make a minimum of bond difference changes
	public int computeInitialLB()
	{
		return new AstarMapping().costEstimate(this);		
	}
	
//	// return the bonds in products, which are equal to bonds in substrates
//	public Collection<Bond> getInducedEdges(Mapping m)
//	{
//		Collection<Atom> mapped_subsatoms = m.getDomain();
//		Collection<Bond> subsbonds = getSubsBonds();
//		Collection<Bond> mapped_subsbonds = new HashSet<Bond>();
//		for (Bond b : subsbonds)
//		{
//			if (mapped_subsatoms.contains(b.getA1()) && mapped_subsatoms.contains(b.getA2()))
//				mapped_subsbonds.add(b);
//		}
//
//		Collection<Bond> mapped_bonds = new HashSet<Bond>();
//		for (Bond b : mapped_subsbonds)
//		{
//			// take the counterparts of this bonds atoms
//			Atom ma1 = m.getImage(b.getA1());
//			Atom ma2 = m.getImage(b.getA2());
//			Bond mb  = ma1.getBond(ma2);
//
//			// if there exists a bond between them (of same type as on left side),
//			// use the same bond
//			if (mb != null && mb.getType() == b.getType())
//				mapped_bonds.add(mb);
//			// otherwise create a new one
//			else
//				mapped_bonds.add( new Bond(ma1, ma2, b.getType(), false));
//		}
//
//		return mapped_bonds;
//	}
//
//	/*
//	 * Finds the neighboring nodes of the set @param big_S in the substrate forest
//	 */
//	public Collection<Atom> getRealNeighbors(Collection<Atom> big_S)
//	{
//		Collection<Bond> s_bonds = getSubsBonds();
//		Collection<Bond> bdrEdges = new HashSet<Bond>();
//		for (Bond b : s_bonds)
//		{
//			if (big_S.contains(b.getA1()) && !big_S.contains(b.getA2()))
//			{
//				bdrEdges.add(b);
//			}
//			if (!big_S.contains(b.getA1()) && big_S.contains(b.getA2()))
//			{
//				bdrEdges.add(b);
//			}
//		}
//		Collection<Atom> border = new HashSet<Atom>();
//		for (Bond b : bdrEdges)
//		{
//			if (!big_S.contains(b.getA1()) && !border.contains(b.getA1()))
//			{
//				border.add(b.getA1());
//			}
//			if (!big_S.contains(b.getA2()) && !border.contains(b.getA2()))
//			{
//				border.add(b.getA2());
//			}
//		}
//
//		return border;
//	}
//
//	
//	/*
//	 * Finds the neighboring nodes of the set @param big_S in the substrate forest
//	 */
//	public Collection<Atom> getNeighbors(Collection<Atom> big_S)
//	{
//		Collection<Atom> s_atoms = getSubsAtoms();
//		Collection<Bond> s_bonds = getSubsBonds();
//		Collection<Bond> bdrEdges = new HashSet<Bond>();
//		for (Bond b : s_bonds)
//		{
//			if (big_S.contains(b.getA1()) && !big_S.contains(b.getA2()))
//			{
//				bdrEdges.add(b);
//			}
//			if (!big_S.contains(b.getA1()) && big_S.contains(b.getA2()))
//			{
//				bdrEdges.add(b);
//			}
//		}
//		Collection<Atom> border = new HashSet<Atom>();
//		for (Bond b : bdrEdges)
//		{
//			if (!big_S.contains(b.getA1()) && !border.contains(b.getA1()))
//			{
//				border.add(b.getA1());
//			}
//			if (!big_S.contains(b.getA2()) && !border.contains(b.getA2()))
//			{
//				border.add(b.getA2());
//			}
//		}
//		if (border.isEmpty())
//			border = Atom.setMinus(s_atoms, big_S);
//		return border;
//	}

	public List<Molecule> getProds()
	{
		return this.products;
	}

	public Collection<Atom> getProdsAtoms()
	{
		if (this.prodatoms.size() == 0)
		{
			for (Molecule m : this.products)
				this.prodatoms.addAll(m.getAtoms());
		}

		return this.prodatoms;
	}

	public Collection<Bond> getProdsBonds()
	{
		if (this.prodbonds.size() == 0)
		{
			for (Molecule m : this.products)
				this.prodbonds.addAll(m.getBonds());
		}

		return this.prodbonds;
	}

	public String getReaction()
	{
		return this.reaction;
	}

	private String[] getReactionString()
	{
		String[] ret = null;
		try
		{
			Scanner sc = new Scanner(new File(GlobalOptions.reacfile));
			while (sc.hasNextLine())
			{
				String str = sc.nextLine();
				ret = str.split(" ");
				if (ret[0].equals(this.id))
					break;
			}
			sc.close();
		} catch (Exception e)
		{
			System.out.println(e.getMessage());
			return null;
		}

		return ret;

	}

	public List<Molecule> getSubs()
	{
		return this.substrates;
	}

	public Collection<Atom> getSubsAtoms()
	{
		if (this.subatoms.size() == 0)
		{
			for (Molecule m : this.substrates)
				this.subatoms.addAll(m.getAtoms());
		}

		return this.subatoms;
	}

	public Collection<Bond> getSubsBonds()
	{
		if (this.subbonds.size() == 0)
		{
			for (Molecule m : this.substrates)
				this.subbonds.addAll(m.getBonds());
		}

		return this.subbonds;
	}

	public void insertMapping(HashMapping m)
	{
		this.mappings.add(m);
	}

	public void printMappingToFile(String path, long time)
	{
		try
		{
			new File(path).mkdir();
			
			BufferedWriter out = new BufferedWriter(new FileWriter(path + this.id + ".txt"));
			int i = 1;
			out.write(printReactionInfo(time));
			out.newLine();

			for (Mapping m : this.mappings)
			{
				if (m == null)
					System.out.println("NULL");
				
				out.write("MAPPING " + (i++));
				out.newLine();
				out.write(m.printMapping());
				
				// we don't use 'BONDS OF REACTION GRAPH' section for anything really...
				out.write("BONDS OF REACTION GRAPH");
				out.newLine();
				out.write("ATOM1\tMOL\tATOM2\tMOL\tORDER\tTYPE");
				out.newLine();
				
				String temp = "";
				String temp2 = "";
				HashSet<String> used = new HashSet<String>();
				
				for (Atom a1 : this.getSubsAtoms())
				{
					for (Atom a2 : this.getSubsAtoms())
					{
						Bond lhs = a1.getBond(a2);
						Bond rhs = m.getImage(a1).getBond( m.getImage(a2) );
						
						int a1row = a1.getMolecularId()+1;
						int a2row = a2.getMolecularId()+1;
						int a1mol = a1.getParent().getMolNum()+1;
						int a2mol = a2.getParent().getMolNum()+1;
						
						if (lhs != null && rhs != null && lhs.getType() == rhs.getType()) // match
						{
							temp = a1row + "\t" + a1mol + "\t" + a2row + "\t" + a2mol + "\t" + lhs.getType() + "\t" + "0" + "\n";
							temp2 = a2row + "\t" + a2mol + "\t" + a1row + "\t" + a1mol + "\t" + lhs.getType() + "\t" + "0" + "\n";
						}
						else if (lhs != null && rhs != null && lhs.getType() != rhs.getType()) // mismatch
						{
							if (GlobalOptions.bc == false)
							{
								temp = a1row + "\t" + a1mol + "\t" + a2row + "\t" + a2mol + "\t" + lhs.getType()+"/"+rhs.getType() + "\t" + "~0" + "\n";
								temp2 = a2row + "\t" + a2mol + "\t" + a1row + "\t" + a1mol + "\t" + lhs.getType()+"/"+rhs.getType() + "\t" + "~0" + "\n";
							}
							else
							{
								temp = a1row + "\t" + a1mol + "\t" + a2row + "\t" + a2mol + "\t" + lhs.getType() + "\t" + "-1" + "\n";
								temp += a1row + "\t" + a1mol + "\t" + a2row + "\t" + a2mol + "\t" + rhs.getType() + "\t" + "1" + "\n";
								temp2 = a2row + "\t" + a2mol + "\t" + a1row + "\t" + a1mol + "\t" + lhs.getType() + "\t" + "-1" + "\n";
								temp2 += a2row + "\t" + a2mol + "\t" + a1row + "\t" + a1mol + "\t" + rhs.getType() + "\t" + "1" + "\n";
							}
						}
						else if (lhs == null && rhs != null) // new bond
						{	
							temp = a1row + "\t" + a1mol + "\t" + a2row + "\t" + a2mol + "\t" + rhs.getType() + "\t" + "1" + "\n";
							temp2 = a2row + "\t" + a2mol + "\t" + a1row + "\t" + a1mol + "\t" + rhs.getType() + "\t" + "1" + "\n";
						}
						else if (lhs != null && rhs == null) // cleaved bond
						{	
							temp = a1row + "\t" + a1mol + "\t" + a2row + "\t" + a2mol + "\t" + lhs.getType() + "\t" + "-1" + "\n";
							temp2 = a2row + "\t" + a2mol + "\t" + a1row + "\t" + a1mol + "\t" + lhs.getType() + "\t" + "-1" + "\n";
						}
						else // lhs == null && rhs == null
							continue;

						
						if (!used.contains(temp) && !used.contains(temp2))
						{
							out.write(temp);
							used.add(temp);
						}
					}
				}
				
				out.newLine();

			}
			out.close();
		}
		catch (IOException e)
		{
			System.out.println("Ioexception! " + e.getMessage());
		}
		
		System.out.println("Reaction file for " + this.id + " written succesfully.");
	}

	public String printReactionInfo(long time)
	{
		StringBuffer sbf = new StringBuffer();
		sbf.append("REACTION " + this.id);
		sbf.append("\n");
		Iterator<Molecule> it = this.substrates.iterator();
		sbf.append("EQUATION ");
		while (it.hasNext())
		{
			sbf.append(it.next().getId() + "\t");
			if (it.hasNext())
				sbf.append(" + ");
		}
		sbf.append("=>\t");
		it = this.products.iterator();
		while (it.hasNext())
		{
			sbf.append(it.next().getId() + "\t");
			if (it.hasNext())
				sbf.append(" + ");
		}
		sbf.append("\n");
		it = this.substrates.iterator();
		sbf.append("INDICES\t ");
		while (it.hasNext())
		{
			sbf.append(it.next().getMolNum()+1 + "\t");
			if (it.hasNext())
				sbf.append(" + ");
		}
		sbf.append("\t=>\t");
		it = this.products.iterator();
		while (it.hasNext())
		{
			sbf.append(it.next().getMolNum()+1 + "\t");
			if (it.hasNext())
				sbf.append(" + ");
		}
		
		// t�h�n v�liin cost
		sbf.append("\nCOST " + this.getFirstMapping().getFCost());
		sbf.append("\nTIME " + (time / 1000.0));
		
		if (GlobalOptions.alg == "mcs")
			sbf.append("\nR 10000");
		
		sbf.append("\n");
		return sbf.toString();
	}

	// reads the reaction based on its 'id'
	// First goes through 'kegg_reactions.txt' to find substrates and products
	//  and then reads them in as Molecule-objects
	public void read() throws Exception
	{
//		Molecule.resetMolCounter();
		
		// get reaction str based on the 'id' from 'kegg_reactions.txt'
		String[] args = getReactionString();

		int begin = 2;
		int j = -1;
		for (int i = begin; i < args.length; i++)
		{
			if (args[i].equals("p"))
				j = i;
		}
		if (j == -1)
			throw new Exception("no 'p' in reactiondescription.");

		int len = args.length;
		List<String> lsubs = new ArrayList<String>();
		List<String> lprods = new ArrayList<String>();
		for (int i = begin; i < j; i++)
			lsubs.add(args[i]);
		for (int i = j + 1; i < len; i++)
			lprods.add(args[i]);

		int molindex = 0;
		int atomindex = 0;
		
		for (String s : lsubs)
		{
			Molecule m = new Molecule(molindex++, s);
			try
			{
				m.read(s + ".mol", GlobalOptions.molpath, atomindex);
			} catch (IOException e)
			{
				throw new Exception(GlobalOptions.molpath + s + ".mol file not found");
			}

			if (!m.getAtoms().isEmpty())
			{
//				for (Atom a : m.getAtoms())
//					a.setReactionalId(atomindex++);
				this.substrates.add(m);
				atomindex += m.getAtoms().size();
			}
			else
				molindex--; // re-use index
		}

		for (String s : lprods)
		{
			Molecule m = new Molecule(molindex++, s);
			try
			{
				m.read(s + ".mol", GlobalOptions.molpath, atomindex);
			} catch (IOException e)
			{
				throw new Exception(GlobalOptions.molpath + s + ".mol file not found");
			}

			if (!m.getAtoms().isEmpty())
			{	
//				for (Atom a : m.getAtoms())
//					a.setReactionalId(atomindex++);
				this.products.add(m);
				atomindex += m.getAtoms().size();
			}
			else
				molindex--;
		}
		
		// Larger molecules first
		Collections.sort(this.substrates);
		Collections.sort(this.products);
		
		this.molsNum = substrates.size() + products.size();
		
		this.reOrder();
		this.computeMaxDistAtoms();
		this.reOrderBFS();
		
		this.getSubsAtoms();
		this.getSubsBonds();
		this.getProdsAtoms();
		this.getProdsBonds();
	}

	public void setId(String id)
	{
		this.id = id;
	}

	public void setProds(List<Molecule> prods)
	{
		this.products = prods;
	}

	public void setReaction(String reac)
	{
		this.reaction = reac;
	}

	public void setSubs(List<Molecule> subs)
	{
		this.substrates = subs;
	}

//	// goes through the mapping and returns bonds which have changed type
//	public Collection<Bond> getBondChanges(Mapping m)
//	{
//		Collection<Bond> result = new HashSet<Bond>();
//
//		// get induced bonds
//		// get product bonds
//		// take differing bonds
//		
//		Collection<Bond> induced = this.getInducedEdges(m);
//		Collection<Bond> prods = this.getProdsBonds();
//		
//		Collection<Bond> intersection = new HashSet<Bond>();
//		
//		intersection.addAll(induced);
//		intersection.retainAll(prods);
//
//		Collection<Bond> set1 = new HashSet<Bond>();
//		set1.addAll(induced);
//		set1.removeAll(prods);
//		
//		return result;
//	}
//	
//	public Collection<Bond> getBondAdditions(Mapping m)
//	{
//		Collection<Bond> result = new HashSet<Bond>();
//		
//		return result;
//	}
//	
//	public Collection<Bond> getBondRemovals(Mapping m)
//	{
//		Collection<Bond> result = new HashSet<Bond>();
//		
//		return result;
//	}
//	
//	public void writeGraphML()
//	{
//		writeGraphML(GlobalOptions.outputdir, this.id);
//	}
//
//	public void writeGraphML(String path, String fn)
//	{
//		String xml = this.toGraphML();
//
//		try
//		{
//			BufferedWriter out = new BufferedWriter(new FileWriter(path + fn + ".graphml"));
//
//			out.write(xml);
//			
//			out.close();
//		} catch (IOException e)
//		{
//			System.out.println(e.getMessage());
//		}
//	}
//	
//	// Creates GraphML representation of the reaction using the mapping
//	public String toGraphML()
//	{
//		StringBuffer sbf = new StringBuffer();
//		sbf.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
//		sbf.append("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n");
//		sbf.append("    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
//		sbf.append("    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n");
//		sbf.append("    http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n");
//		sbf.append("<graph id=\"G\" edgedefault=\"undirected\">\n");
//
//		sbf.append(toGraphML(this.mappings.get(0).getRange(),
//				getInducedEdges(this.mappings.get(0)), 1, "M"));
//		sbf.append(toGraphML(getProdsAtoms(), getProdsBonds(), 2, "P"));
//		sbf.append("</graph>\n");
//		sbf.append("</graphml>\n");
//
//		return sbf.toString();
//	}
//
//	private String toGraphML(Collection<Atom> atoms, Collection<Bond> bonds, int clus,
//			String side)
//	{
//		StringBuffer sbf = new StringBuffer();
//
//		for (Atom a : atoms)
//		{
//			sbf.append("<node id=\"");
//			sbf.append("" + side + a.getId());
//			sbf.append("\">\n");
//			sbf.append("\t<data key=\"text\">");
//			sbf.append("" + a.getSymbol() + a.getId() + " (" + side + ")");
//			sbf.append("</data>\n");
//			sbf.append("\t<data key=\"clusterID\">");
//			sbf.append(clus);
//			sbf.append("</data>\n");
//			sbf.append("\t<data key=\"shape\">");
//			sbf.append("Ellipse");
//			sbf.append("</data>\n");
//			sbf.append("\t<data key=\"color\">");
//			if (a.getSymbol().equals("C"))
//				sbf.append("199 199 199");
//			else if (a.getSymbol().equals("O"))
//				sbf.append("255 0 0");
//			else if (a.getSymbol().equals("N"))
//				sbf.append("0 255 0");
//			else if (a.getSymbol().equals("S"))
//				sbf.append("255 255 0");
//			else if (a.getSymbol().equals("P"))
//				sbf.append("0 0 255");
//			else
//				sbf.append("255 255 255");
//			sbf.append("</data>\n");
//			sbf.append("</node>\n");
//		}
//		for (Bond b : bonds)
//		{
//			sbf.append("<edge source=\"");
//			sbf.append("" + side + b.getA1().getId());
//			sbf.append("\" target=\"");
//			sbf.append("" + side + b.getA2().getId());
//			sbf.append("\">\n");
//
//			sbf.append("\t<data key=\"width\">");
//			sbf.append(1);
//			sbf.append("</data>");
//
//			sbf.append("</edge>\n");
//
//		}
//		return sbf.toString();
//	}
	
	public void computeMaxDistAtoms()
	{
		for (Molecule m : this.getSubs())
			m.computeMaxDistAtom();
			
		for (Molecule m : this.getProds())
			m.computeMaxDistAtom();
	}
	
	// Arbitrary atom ordering
	public void reOrder()
	{
		int i = 0;
		for (Atom a : getSubsAtoms())
			a.setOrder(i++);

		i = 0;
		for (Atom a : getProdsAtoms())
			a.setOrder(i++);
	}
	
	public void resetOrdering(Collection<Molecule> mols)
	{
		for (Molecule m : mols)
			for (Atom a : m.getAtoms())
				a.setOrder(-1);
	}
	
	public void reOrderBFS()
	{
		reOrderBFS(this.getSubs());
		reOrderBFS(this.getProds());
	}
	
	public void reOrderDFS()
	{
		reOrderDFS(this.getSubs());
		reOrderDFS(this.getProds());
	}
	
	public void reOrderBFS(Collection<Molecule> mols)
	{
		resetOrdering(mols);
		
		LinkedList<Atom> Q = new LinkedList<Atom>();
		int i = 0;
		Atom node;
	
		for (Molecule m : mols)
		{
			Q.addFirst(m.getMaxDistAtom());
				
			while (!Q.isEmpty())
			{
				node = Q.removeLast();
				
				if (node.getOrder() >= 0)
					continue;
				
				node.setOrder(i++);
				
				for (Atom ne : node.getAtomNeighbors())
					if (ne.getOrder() < 0)
						Q.addFirst(ne);
			}
		}
		
		for (Molecule m : mols)
			for (Atom a : m.getAtoms())
				if (a.getOrder() < 0)
					a.setOrder(i++);
	}

	public void reOrderDFS(Collection<Molecule> mols)
	{
		resetOrdering(mols);

		LinkedList<Atom> Q = new LinkedList<Atom>();
		int i = 0;
		Atom node;
		Atom start = mols.iterator().next().getMaxDistAtom();
	
		Q.addFirst(start);
		
		for (Molecule m : mols)
		{
			if (!m.getAtoms().contains(start))
				Q.add(m.getMaxDistAtom());
				
			while (!Q.isEmpty())
			{
				node = Q.removeFirst();
				
				if (node.getOrder() >= 0)
					continue;
				
				node.setOrder(i++);
				
				for (Atom ne : node.getAtomNeighbors())
					if (ne.getOrder() < 0)
						Q.addFirst(ne);
			}
		}
		
		for (Molecule m : mols)
			for (Atom a : m.getAtoms())
				if (a.getOrder() < 0)
					a.setOrder(i++);
	}
	
	
	public int traverse(Atom start)
	{
		int[] dist = new int[start.getParent().getAtoms().size()];
		LinkedList<Atom> Q = new LinkedList<Atom>();
		Q.addFirst(start);
		Atom node;
		dist[start.getOrder()] = 0;
		
		while (!Q.isEmpty())
		{
			node = Q.removeLast();
			
			for (Atom ne : node.getAtomNeighbors())
			{
				if (dist[ne.getOrder()] == 0)
				{
					Q.addFirst(ne);
					dist[ne.getOrder()] = dist[node.getOrder()] + 1;
				}
			}
		}
		
		int max = 0;
		for (int x : dist)
			if (x > max)
				max = x;
		
		return max;
	}


	public void constructFeatureTables() throws Exception
	{
		int submaxid = 0;
		int prodmaxid = 0;
		int maxrow = 0;
		
		for (Atom a : this.getSubsAtoms())
		{
			if (a.getOrder() > submaxid)
				submaxid = a.getOrder();
		}
		for (Atom a : this.getProdsAtoms())
		{			
			if (a.getOrder() > prodmaxid)
				prodmaxid = a.getOrder();
		}

		submaxid++;
		prodmaxid++;
		
		// TODO should we add +1 to size()?
		for (Molecule m : this.getSubs())
			if (m.getAtoms().size() > maxrow)
				maxrow = m.getAtoms().size();
		for (Molecule m : this.getProds())
			if (m.getAtoms().size() > maxrow)
				maxrow = m.getAtoms().size();
		
		
	
		
//		wiener = new double[submaxid][prodmaxid];
		
		// kahden atomin v�linen piirrevektori on konteksteittain eri suuntien samuus
		// atom 3 <-> atom 9': ( (1,1),(1,1),(0,0),(0,0) )
		//                   => (1,1,0,0)
		
		morgan = new double[11][submaxid][prodmaxid];
		wiener = new double[11][submaxid][prodmaxid];
		ad = new double[11][submaxid][prodmaxid];
		bd = new double[11][submaxid][prodmaxid];
//		ring = new double[7][submaxid][prodmaxid];
		
		featcost = new double[submaxid][prodmaxid];
		similarity = new double[submaxid][prodmaxid];
		
		int[][][][] morganfeat = new int[molsNum][11][maxrow][5];
		int[][][][] wienerfeat = new int[molsNum][11][maxrow][5];
		String[][][][] adfeat = new String[molsNum][11][maxrow][5];
		String[][][][] bdfeat = new String[molsNum][11][maxrow][5];
		
		
		// luetaan feat-linet mol-taulukkoihin
		// molecule -> context -> atomid -> value
		
		for (Molecule m : this.getSubs())
		{
			String featfile = GlobalOptions.featpath + m.getId() + ".feat";
			String line;
			Scanner sc = null;
			
			try {
				sc = new Scanner(new File(featfile));
			} catch (Exception e)
			{
				System.out.println("File not found: " + featfile);
//				GlobalOptions.plusplus = false;
				throw e;
			}
		
			while (sc.hasNextLine())
			{
//				C00040:13:BOND_DISTRIBUTION:1=(111111222,7=112,17=12,18=1)	
//				C00040:2:BOND_DISTRIBUTION:4=(111111111111111111111111111111111111122222222,1=11111111111112222222,5=11111111111111111,6=11111111111111111112)
//				C00040:12:ATOM_DISTRIBUTION:2=(CCCO,6=CCC)

				line = sc.nextLine();
				
				String[] words = line.split(":");
				
				String feat = words[2];
				
				if (feat.equals("RING"))
					continue;
				
				int context = Integer.parseInt(words[3].split("=")[0]);
				int atomnum = Integer.parseInt(words[1]);
//				String value = words[3].split("\\(")[1].split("\\)")[0].split(",")[0];

				String[] value;

				if (words[3].contains("()"))
				{
					value = new String[1];
					value[0] = "";
				}
				else
					value = words[3].split("\\(")[1].split("\\)")[0].split(",[0-9]+=");
				
				
				// huge kludge !!!
				// sometimes feat-line is "... WIENER, 1=(0,)", in which case
				// value[0] is "0,"
				if (value[0].contains(","))
					value[0] = value[0].substring(0,value[0].length()-1);
				
				Arrays.sort(value, 1, value.length);
				
				
				if (words[2].equals("WIENER"))
				{
					int[] intvalue = new int[value.length];
					int i = 0;
					for (String s : value)
						intvalue[i++] = Integer.parseInt(s);

//					a.feats.set("wiener", context, value)
					
					wienerfeat[m.getMolNum()][context][atomnum-1] = intvalue;
				}
				else if (words[2].equals("MORGAN"))
				{
					int[] intvalue = new int[value.length];
					int i = 0;
					for (String s : value)
						intvalue[i++] = Integer.parseInt(s);
					
					morganfeat[m.getMolNum()][context][atomnum-1] = intvalue;
				}
				else if (words[2].equals("ATOM_DISTRIBUTION"))
				{
					adfeat[m.getMolNum()][context][atomnum-1] = value;
				}
				else if (words[2].equals("BOND_DISTRIBUTION"))
				{
					bdfeat[m.getMolNum()][context][atomnum-1] = value;
				}
			}
		}
		for (Molecule m : this.getProds())
		{
			String featfile = GlobalOptions.featpath + m.getId() + ".feat";
			String line;
			Scanner sc = null;
			
			try {
				sc = new Scanner(new File(featfile));
			} catch (Exception e)
			{
				System.out.println("File not found: " + featfile);
//				GlobalOptions.plusplus = false;
				throw e;
			}
		
			while (sc.hasNextLine())
			{
//				C00040:13:BOND_DISTRIBUTION:1=(111111222,7=112,17=12,18=1)	
//				C00040:2:BOND_DISTRIBUTION:4=(111111111111111111111111111111111111122222222,1=11111111111112222222,5=11111111111111111,6=11111111111111111112)
//				C00040:12:ATOM_DISTRIBUTION:2=(CCCO,6=CCC)

				line = sc.nextLine();
				
				String[] words = line.split(":");
				
				String feat = words[2];
				
				if (feat.equals("RING"))
					continue;
				
				int context = Integer.parseInt(words[3].split("=")[0]);
				int atomnum = Integer.parseInt(words[1]);
//				String value = words[3].split("\\(")[1].split("\\)")[0].split(",")[0];
				String[] value;

				if (words[3].contains("()"))
				{
					value = new String[1];
					value[0] = "";
				}
				else
					value = words[3].split("\\(")[1].split("\\)")[0].split(",[0-9]+=");
				
				// huge kludge !!!
				// sometimes feat-line is "... WIENER, 1=(0,)", in which case
				// value[0] is "0,"
				if (value[0].contains(","))
					value[0] = value[0].substring(0,value[0].length()-1);
				
				Arrays.sort(value, 1, value.length);

				if (words[2].equals("WIENER"))
				{
					int[] intvalue = new int[value.length];
					int i = 0;
					for (String s : value)
						intvalue[i++] = Integer.parseInt(s);

					wienerfeat[m.getMolNum()][context][atomnum-1] = intvalue;
				}
				else if (words[2].equals("MORGAN"))
				{
					int[] intvalue = new int[value.length];
					int i = 0;
					for (String s : value)
						intvalue[i++] = Integer.parseInt(s);
					
					morganfeat[m.getMolNum()][context][atomnum-1] = intvalue;
				}
				else if (words[2].equals("ATOM_DISTRIBUTION"))
				{
					adfeat[m.getMolNum()][context][atomnum-1] = value;
				}
				else if (words[2].equals("BOND_DISTRIBUTION"))
				{
					bdfeat[m.getMolNum()][context][atomnum-1] = value;
				}
			}
		}
		
		// molecule -> context -> atomid -> value
		
		for (Atom a1 : this.getSubsAtoms())
		{
			for (Atom a2 : this.getProdsAtoms())
			{
				for (int i = 1; i <= 10; i++)
				{
					int[] a1wiener = wienerfeat[a1.getParent().getMolNum()][i][a1.getMolecularId()];
					int[] a2wiener = wienerfeat[a2.getParent().getMolNum()][i][a2.getMolecularId()];

					if (a1wiener[0] ==  0 || a2wiener[0] == 0)
					{
						wiener[i][a1.getOrder()][a2.getOrder()] = 1.0;
					}
					else
					{
						if (a1wiener[0] == a2wiener[0])
							wiener[i][a1.getOrder()][a2.getOrder()] = 0.0;
						else
						{
							// TODO t�h�n kohtaan k�yt�v� l�pi a1 ja a2 wiener-taulukot ja m�ts�tt�v� kuinka
							// suuri osuus eri suunnista m�ts��
							
							// loopattava temp1[1:end] ja temp2[1:end]
							int match_branches = 0;
							double branches = Math.max(a1wiener.length, a2wiener.length) - 1;
							int a = 1, b = 1;
							while (a < a1wiener.length && b < a2wiener.length)
							{
								if (a1wiener[a] == a2wiener[b])
								{
									match_branches++;
									a++;
									b++;
								}
								else if (a1wiener[a] < a2wiener[b])
									a++;
								else
									b++;
							}
							
							wiener[i][a1.getOrder()][a2.getOrder()] = 1.0 - (match_branches / branches);
						}
					}
					
					
					String[] a1ad = adfeat[a1.getParent().getMolNum()][i][a1.getMolecularId()];
					String[] a2ad = adfeat[a2.getParent().getMolNum()][i][a2.getMolecularId()];

					if (a1ad[0] == null || a2ad[0] == null)
					{
						ad[i][a1.getOrder()][a2.getOrder()] = 1.0;
					}
					else
					{
						if (a1ad[0].equals(a2ad[0]))
							ad[i][a1.getOrder()][a2.getOrder()] = 0.0;
						else
						{
							// TODO t�h�n kohtaan k�yt�v� l�pi a1 ja a2 wiener-taulukot ja m�ts�tt�v� kuinka
							// suuri osuus eri suunnista m�ts��
							
							// loopattava temp1[1:end] ja temp2[1:end]
							int match_branches = 0;
							double branches = Math.max(a1ad.length, a2ad.length) - 1;
							int a = 1, b = 1;
							while (a < a1ad.length && b < a2ad.length)
							{
								if (a1ad[a].equals(a2ad[b]))
								{
									match_branches++;
									a++;
									b++;
								}
								else if (a1ad[a].compareTo(a2ad[b]) < 0)
									a++;
								else
									b++;
							}
							
							ad[i][a1.getOrder()][a2.getOrder()] = 1.0 - (match_branches / branches);
						}
					}

					int[] a1morgan = morganfeat[a1.getParent().getMolNum()][i][a1.getMolecularId()];
					int[] a2morgan = morganfeat[a2.getParent().getMolNum()][i][a2.getMolecularId()];

					if (a1morgan[0] ==  0 || a2morgan[0] == 0)
					{
						morgan[i][a1.getOrder()][a2.getOrder()] = 1.0;
					}
					else
					{
						if (a1morgan[0] == a2morgan[0])
							morgan[i][a1.getOrder()][a2.getOrder()] = 0.0;
						else
						{
							// TODO t�h�n kohtaan k�yt�v� l�pi a1 ja a2 wiener-taulukot ja m�ts�tt�v� kuinka
							// suuri osuus eri suunnista m�ts��
							
							// loopattava temp1[1:end] ja temp2[1:end]
							int match_branches = 0;
							double branches = Math.max(a1morgan.length, a2morgan.length) - 1;
							int a = 1, b = 1;
							while (a < a1morgan.length && b < a2morgan.length)
							{
								if (a1morgan[a] == a2morgan[b])
								{
									match_branches++;
									a++;
									b++;
								}
								else if (a1morgan[a] < a2morgan[b])
									a++;
								else
									b++;
							}
							
							morgan[i][a1.getOrder()][a2.getOrder()] = 1.0 - (match_branches / branches);
						}
					}
					
					
//					wiener[i][a1.getOrder()][a2.getOrder()] = (wienerfeat[a1.getParent().getMolNum()-1][i][a1.getRownumber()-1] == wienerfeat[a2.getParent().getMolNum()-1][i][a2.getRownumber()-1]) ? 0 : 1;
//					ad[i][a1.getOrder()][a2.getOrder()] = (adfeat[a1.getParent().getMolNum()-1][i][a1.getRownumber()-1].equals( adfeat[a2.getParent().getMolNum()-1][i][a2.getRownumber()-1] ) ) ? 0 : 1;
//					bd[i][a1.getOrder()][a2.getOrder()] = (bdfeat[a1.getParent().getMolNum()-1][i][a1.getRownumber()-1].equals( bdfeat[a2.getParent().getMolNum()-1][i][a2.getRownumber()-1] ) ) ? 0 : 1;
				}
			}
		}
		
		for (Atom a1 : getSubsAtoms())
		{
			for (Atom a2 : getProdsAtoms())
			{
				double c = 0.0;
				int a1n = a1.getOrder();
				int a2n = a2.getOrder();
				
				int similar_count = 0;
				
				for (int k = 1; k <= 10; k++)
				{
					similar_count  += (ad[k][a1n][a2n] == 0.0) ? 1 : 0;
					similar_count  += (wiener[k][a1n][a2n] == 0.0) ? 1 : 0;
					similar_count  += (morgan[k][a1n][a2n] == 0.0) ? 1 : 0;
					
					similarity[a1n][a2n] += (ad[k][a1n][a2n] == 0.0) ? 1.0 : 0.0;
					similarity[a1n][a2n] += (wiener[k][a1n][a2n] == 0.0) ? 1.0 : 0.0;
					similarity[a1n][a2n] += (morgan[k][a1n][a2n] == 0.0) ? 1.0 : 0.0;
				
					c += wiener[k][a1.getOrder()][a2.getOrder()];
					c += morgan[k][a1.getOrder()][a2.getOrder()];
					c += ad[k][a1.getOrder()][a2.getOrder()];
					
					if (GlobalOptions.bc == true)
						c += bd[k][a1.getOrder()][a2.getOrder()];
				}

				similarity[a1n][a2n] /= 40.0;
				
				similarity[a1n][a2n] = similar_count / 40.0;
				
				double total = (GlobalOptions.bc == false) ? 30.0 : 40.0;
				
				featcost[a1.getOrder()][a2.getOrder()] = c / total;
				
			}
		}
	}
	
	/*
	 * l�ht�kohtaisesti halutaan [0..1] v�lilt� featcost f(a,b)
	 * t�m� muodostuu ad,wiener,bd featureista tasoilla 1,2,3,4 a:n ja b:n v�lill�.
	 * 
	 * lis�n� entrooppinen h�ss�kk� joka ottaa huomioon eri suunnat. eli 
	 */
	
	
	public double normalizedFeatCost(Atom a1, Atom a2)
	{
		return featcost[a1.getOrder()][a2.getOrder()];
	}
	
	public double featCost(Atom a1, Atom a2)
	{
		double cost = 0.0;
		
		for (int k = 1; k <= 10; k++)
		{
			cost += 0.0078125 * wiener[k][a1.getOrder()][a2.getOrder()];
			cost += 0.0078125 * morgan[k][a1.getOrder()][a2.getOrder()];
			cost += 0.0078125 * ad[k][a1.getOrder()][a2.getOrder()];
		
			if (GlobalOptions.bc == true)
				cost += 0.0078125 * bd[k][a1.getOrder()][a2.getOrder()];
		}
		
		return cost;
	}
	
	public boolean featDiff(Atom a1, Atom a2, int context)
	{
		return wiener[context][a1.getOrder()][a2.getOrder()] == 1;
	}
	
	class Distribution
	{
		private Collection<Double> dist = null;
		
		public Distribution(double[] dist)
		{
			this.dist = new ArrayList<Double>();
			
			for (int i = 0; i < dist.length; i++)
				this.dist.add(dist[i]);
			
			normalize();
		}
		
		public Distribution(Collection<Double> dist)
		{
			this.dist = dist;
			
			normalize();
		}
		
		// Normalize distribution to sum to 1
		private void normalize()
		{
			double sum = 0.0;
			
			for (Double d : this.dist)
				sum += d;
			
			for (Double d : this.dist)
				d /= sum;
		}
		
		public double entropy()
		{
			double e = 0.0;
			
			for (Double d : this.dist)
				e -= (d < 0.001) ? 0.0 : Math.log10(d) * d;
			
			return e;
		}
	}
	
	
	// remove atoms from reaction
	public void remove(Mapping map)
	{
		// remove from subbonds and subatoms
		// remove also from molecule.bonds and molecule.atoms
		// remove also from atom its neighbors
		for (Atom a : map.getDomain())
		{
			removedAtoms.add(a);
			
			for (Atom ne : a.getAtomNeighbors())
			{
				removedBonds.add( a.getBond(ne) );
				ne.removeNeighbor(a);
			}
			
			this.subbonds.removeAll(a.getBondNeighbors());
			this.subatoms.remove(a);

			for (Molecule s : this.substrates)
			{
				s.getBonds().removeAll(a.getBondNeighbors());
				s.getAtoms().remove(a);
			}
		}
		
		for (Atom a : map.getRange())
		{
			removedAtoms.add(a);

			for (Atom ne : a.getAtomNeighbors())
			{
				removedBonds.add( a.getBond(ne) );
				ne.removeNeighbor(a);
			}
			
			this.prodbonds.removeAll(a.getBondNeighbors());
			this.prodatoms.remove(a);

			for (Molecule p : this.products)
			{				
				p.getBonds().removeAll(a.getBondNeighbors());
				p.getAtoms().remove(a);
			}
		}
	}
	
	public void restore()
	{
		for (Atom a : removedAtoms)
		{
			a.getParent().addAtom(a);
			
			if (substrates.contains(a.getParent()))
				subatoms.add(a);
			else
				prodatoms.add(a);
		}
		
		for (Bond b : removedBonds)
		{
			b.getParent().addBond(b);
			b.getA1().addNeighbor(b, b.getA2());
			b.getA2().addNeighbor(b, b.getA1());
			
			if (substrates.contains(b.getParent()))
				subbonds.add(b);
			else
				prodbonds.add(b);
		}
		
		removedAtoms.clear();
		removedBonds.clear();
	}
	
	public void computeDistances()
	{
		// compute all pairwise distances of atoms
		distances.clear();
		
		for (Atom a : this.getSubsAtoms())
			distancesFrom(a);
		for (Atom a : this.getProdsAtoms())
			distancesFrom(a);
	}

	// find all 
	private void distancesFrom(Atom a)
	{
		LinkedList<Atom> Q = new LinkedList<Atom>();
		Q.addFirst(a);
		Atom node;
		
		distances.put(a, new HashMap<Atom, Integer>());
		distances.get(a).put(a, 0);
		
		while (!Q.isEmpty())
		{
			node = Q.removeLast();
			
			for (Atom ne : node.getAtomNeighbors())
			{
				if (!distances.get(a).containsKey(ne))
				{
					Q.addFirst(ne);
					distances.get(a).put(ne, distances.get(a).get(node) + 1 );
				}
			}
		}
	}
	
	public int distance(Atom a, Atom b)
	{
		Map<Atom, Integer> dist = distances.get(a);
		
		if (dist.containsKey(b))
			return dist.get(b);
		else
			return -1;
	}

	public void writeProgressInfo(int[] lbhistory, int[] ubhistory)
	{
		try
		{
			BufferedWriter out = new BufferedWriter(new FileWriter(GlobalOptions.outputdir + this.id + ".progress"));
			
			out.write("REACTION " + this.id + "\n");
			
			for (int i = 0; i < ubhistory.length; i++)
				out.write(i + "\t" + lbhistory[i] + "\t" + ubhistory[i] + "\n");
			
			out.close();
		} catch (IOException e)
		{
			System.out.println("Ioexception! " + e.getMessage());
		}
		
		System.out.println("Progress information written to " + GlobalOptions.outputdir + id + ".progress");
	}

//	// return the bonds of 'bonds' which are not between two atoms of 'atoms'
//	public static Collection<Bond> getResidualEdges(Collection<Bond> bonds,
//			Collection<Atom> set)
//	{
//		Collection<Bond> resbonds = new ArrayList<Bond>();
//		for (Bond b : bonds)
//		{
//			if (!(set.contains(b.getA1()) && set.contains(b.getA2())))
//				resbonds.add(b);
//		}
//	
//		return resbonds;
//	}
//	

	// each atom has a neighbor-string (e.g. "C|COO") representing itself and neighbors
	// compute these strings and the indices attached to them, (effectively hash values)
	public void computeAtomStringIndex()
	{
		Map<String,Integer> neighstrs = new HashMap<String,Integer>();
		int counter = 0;
		
		for (Atom a : subatoms)
		{
			String nstr = a.getNeighString();
			
			if (neighstrs.containsKey(nstr))
				a.setNeighStringInt( neighstrs.get(nstr) );
			else
			{
				a.setNeighStringInt( counter );
				neighstrs.put(nstr, counter);
				counter++;
			}
		}
		
		for (Atom a : prodatoms)
		{
			String nstr = a.getNeighString();
			
			if (neighstrs.containsKey(nstr))
				a.setNeighStringInt( neighstrs.get(nstr) );
			else
			{
				a.setNeighStringInt( counter );
				neighstrs.put(nstr, counter);
				counter++;
			}
		}
		
		Astar.atomstrsize = counter;
	}
	
	// each bond is represented by a string, e.g. "CC" for C-C
	// make an index for these
	public void computeBondStringIndex()
	{
		Map<String,Integer> neighstrs = new HashMap<String,Integer>();
		int counter = 0;
		
		for (Bond b : subbonds)
		{
			String nstr = b.getBondString();
			
			if (neighstrs.containsKey(nstr))
				b.setBondStringIndex( neighstrs.get(nstr) );
			else
			{
				b.setBondStringIndex( counter );
				neighstrs.put(nstr, counter);
				counter++;
			}
		}
		
		for (Bond b : prodbonds)
		{
			String nstr = b.getBondString();
			
			if (neighstrs.containsKey(nstr))
				b.setBondStringIndex( neighstrs.get(nstr) );
			else
			{
				b.setBondStringIndex( counter );
				neighstrs.put(nstr, counter);
				counter++;
			}
		}
		
		Astar.bondstrsize = counter;
	}

	
	// writes all reaction graphs as mol files in 'path'
	public void writeReactionGraphs(String path)
	{
		int i = 0;
		for (Mapping m : this.mappings)
		{
			ReactionGraph rg = new ReactionGraph(this, m);
			rg.writeAsMol(path + this.id + "_" + i + "_f.mol", Direction.FORWARD);
			rg.writeAsMol(path + this.id + "_" + i + "_b.mol", Direction.BACKWARD);
			i++;
		}
		
		System.out.println(i*2 + " reaction graphs for " + this.id + " written succesfully.");
	}
		
		
//		
////		System.out.println("startrgraph");
//		int i = 0;
//		for (Mapping m : this.mappings)
//		{
//			writeReactionGraphMol(m, path + this.id + "_" + i + "_f.mol", +1);
////			System.out.println("f");
//			writeReactionGraphMol(m, path + this.id + "_" + i + "_b.mol", -1);
////			System.out.println("b");
//			i++;
//		}
//		System.out.println("donergraph");
//		System.out.println(i*2 + " reaction graphs for " + this.id + " written succesfully.");
//	}
//	
//	// writes reaction graph as a mol file for both directions
//	private void writeReactionGraphMol(Mapping m, String filename, int dir)
//	{
//		BufferedWriter out;
//
//		// atoms
//		int atomcount = this.getSubsAtoms().size();
//		
////		System.out.println("atomcount " + atomcount);
//		// save atoms in id-order
//		Atom[] atoms = new Atom[atomcount];
//		for (Atom a : this.getSubsAtoms())
//		{
////			System.out.println(a.getReactionalId());
//			atoms[a.getReactionalId()] = a;
//		}
////		System.out.println("x");
//
//		// count number of reaction graph bonds
//		int bondcount = 0;
//		for (Atom a1 : this.getSubsAtoms())
//		{
//			for (Atom a2 : this.getSubsAtoms())
//			{
//				if (a1.getReactionalId() >= a2.getReactionalId())
//					continue;
//				
//				Bond lhs = a1.getBond(a2);
//				Bond rhs = m.getImage(a1).getBond( m.getImage(a2) );
//				
//				if (lhs != null || rhs != null)
//					bondcount++;
//			}
//		}
//		
////		System.out.println("x");
//		String atomstr = Mapper2000.rjust(Integer.toString(atomcount),3," ");
//		String bondstr = Mapper2000.rjust(Integer.toString(bondcount),3," ");
//		
//		try
//		{
//			out = new BufferedWriter(new FileWriter(filename));
//			out.write(this.id + "\n\n\n");
//			out.write(atomstr + bondstr + "  0  0  0  0            999 V2000\n");
//			
//			// atoms
//			for (int k = 0; k < atoms.length; k++)
//				out.write("    0.0000    0.0000    0.0000 " + atoms[k].getSymbol() + "   0  0  0  0  0  0  0  0  0  0  0  0\n");
//			
//			// bonds
//			for (Atom a1 : this.getSubsAtoms())
//			{
//				for (Atom a2 : this.getSubsAtoms())
//				{
//					if (a1.getReactionalId() >= a2.getReactionalId())
//						continue;
//
//					Bond lhs = a1.getBond(a2);
//					Bond rhs = m.getImage(a1).getBond( m.getImage(a2) );
//					
//					if (lhs != null || rhs != null)
//					{
//						out.write( Mapper2000.rjust(Integer.toString(a1.getReactionalId()+1), 3, " ") ); // source id
//						out.write( Mapper2000.rjust(Integer.toString(a2.getReactionalId()+1), 3, " ") ); // target id
//						if (lhs != null && rhs == null) // cleaved bond
//						{
//							if (dir == 1)
//							{
//								out.write(" -1");
//								out.write( Mapper2000.rjust(Integer.toString(lhs.getType()), 3, " ") ); // lhs side bond type
//								out.write("  0");
//							}
//							else
//							{
//								out.write("  1");
//								out.write("  0");
//								out.write( Mapper2000.rjust(Integer.toString(lhs.getType()), 3, " ") ); // lhs side bond type
//							}
//						}
//						else if (lhs == null && rhs != null) // new bond
//						{
//							if (dir == 1)
//							{
//								out.write("  1");
//								out.write("  0");
//								out.write( Mapper2000.rjust(Integer.toString(rhs.getType()), 3, " ") ); // rhs side bond type
//							}
//							else
//							{
//								out.write(" -1");
//								out.write( Mapper2000.rjust(Integer.toString(rhs.getType()), 3, " ") ); // rhs side bond type
//								out.write("  0");
//							}
//						}
//						else if (lhs != null && rhs != null) // no-change
//						{	
//							out.write("  0");
//							out.write( Mapper2000.rjust(Integer.toString(lhs.getType()), 3, " ") ); // lhs side bond type
//							out.write( Mapper2000.rjust(Integer.toString(rhs.getType()), 3, " ") ); // rhs side bond type
//						}
//						out.write( "  0  0\n");
//					}
//				}
//			}
//			
//			out.write("M  END\n");
//			out.close();
//		}
//		catch (IOException e)
//		{
//			System.out.println("X");
//			System.out.println("ioerror: " + e.getMessage());
//			return;
//		}
//	}

	
	// TODO
	// test whether isomorphism algorithm (VF2) works
	// test that reaction graphs work
	// code dir-1/+1 output methods for RG's
	
	
	// tests all mappings for isomorphisms and retains only truely unique mappings
	public void contractIsomorphic()
	{
		List<Mapping> uniquemaps = new ArrayList<Mapping>();
		uniquemaps.add(mappings.get(0)); // prime with first mapping
		
		int i = 0;
		
		for (Mapping m : mappings)
		{
			ReactionGraph rg1 = new ReactionGraph(this, m);
			boolean isomorphic = false;
			
			for (Mapping um : uniquemaps)
			{
				ReactionGraph rg2 = new ReactionGraph(this, um);
				
				Isomorphism iso = new Isomorphism(rg1, rg2);
				if (iso.VF2()) // isomorphic mapping
				{
					isomorphic = true;
					break;
				}
			}
			
			// we've found a new unique mapping
			if (isomorphic == false)
			{
				uniquemaps.add(m);
			}
			
			if (i % 1000 == 0)
			{
				System.out.print(i + "/" + mappings.size() + "(" + uniquemaps.size() + ") ");
				System.out.flush();
			}
			if (i % 10000 == 0)
				System.out.println();
			
				
			i++;
		}
		
		mappings = uniquemaps; // remove non-unique maps
	}

	// product a 2d array of prodatoms:
	// outer array contains bins of atoms (ordered by 'order')
	// inner arrays contain atoms of a single bin (isomorphic)
	// inner arrays are jagged
	public Atom[][] getProdAtomBins()
	{
		Atom[][] arr;
		
		// put atoms by 'order' into singleton bins
		List<List<Atom>> bins = new ArrayList<List<Atom>>();
		for (Atom a : prodatoms)
			bins.add(new ArrayList<Atom>());
		for (Atom a : prodatoms)
			bins.get(a.getOrder()).add(a);
		
		// find sets of atoms with degree 1, shared neighbor and equal atom type
//		List<List<Atom>> groups = new LinkedList<List<Atom>>();
		for (Atom a : prodatoms)
		{
			Map<String,List<Atom>> m = new HashMap<String,List<Atom>>();
			for (Atom ne : a.getAtomNeighbors())
			{
				if (ne.getDegree() == 1)
				{
					String s = ne.getSymbol();
					if (m.containsKey(s))
					{
						m.get(s).add(ne);
					}
					else
					{
						m.put(s, new ArrayList<Atom>());
						m.get(s).add(ne);
					}
				}
			}
			
			// if grouping found, transfer all atoms in the group to the first atoms bin
			for (String s : m.keySet())
			{
				if (m.get(s).size() > 1) // bin found!
				{
					List<Atom> temp = m.get(s);
					
					// add atoms 1..n into bin of 0-atom
					for (int i = 1; i < temp.size(); i++)
						bins.get(temp.get(0).getOrder()).add( temp.get(i) );
						
					// empty bins of 1..n
					for (int i = 1; i < temp.size(); i++)
						bins.get(temp.get(i).getOrder()).clear();
				}
			}
		}
		
		// check the number of non-empty bins
		int bincount = 0;
		for (List<Atom> l : bins)
		{
			if (l.size() > 0)
				bincount++;
		}
		
		// create appropriate size array and put atoms into it
		arr = new Atom[bincount][];
		int i = 0;
		for (List<Atom> l : bins)
		{
			int s = l.size();
			if (s > 0 )
			{
				arr[i] = new Atom[s];
				for (int k = 0; k < s; k++)
				{	
					arr[i][k] = l.get(k);
					l.get(k).bin_id = i;
					l.get(k).bin_offset = k;
				}
				
				i++;
			}
		}
		
		return arr;
	}
}


